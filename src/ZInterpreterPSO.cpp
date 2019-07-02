/* Copyright (C) 2014-2017 Carl Leonardsson
 * Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2019 Viktor Toman
 *
 * This file is part of Nidhugg.
 *
 * Nidhugg is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Nidhugg is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <llvm/IR/Instructions.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/ExecutionEngine/GenericValue.h>
#if defined(HAVE_LLVM_INSTVISITOR_H)
#include <llvm/InstVisitor.h>
#elif defined(HAVE_LLVM_IR_INSTVISITOR_H)
#include <llvm/IR/InstVisitor.h>
#elif defined(HAVE_LLVM_SUPPORT_INSTVISITOR_H)
#include <llvm/Support/InstVisitor.h>
#endif
#if defined(HAVE_LLVM_IR_LLVMCONTEXT_H)
#include <llvm/IR/LLVMContext.h>
#elif defined(HAVE_LLVM_LLVMCONTEXT_H)
#include <llvm/LLVMContext.h>
#endif

#include "ZInterpreterPSO.h"
#include "ZBuilderPSO.h"


static void SetValue(llvm::Value *V, llvm::GenericValue Val, llvm::ExecutionContext &SF) {
  SF.Values[V] = Val;
}


ZInterpreterPSO::ZInterpreterPSO(llvm::Module *M, ZBuilderPSO &TB,
                                 const Configuration &conf)
  : PSOInterpreter(M,TB,conf), TB(TB) {}


ZInterpreterPSO::~ZInterpreterPSO() {}


llvm::ExecutionEngine *ZInterpreterPSO::create(llvm::Module *M, ZBuilderPSO &TB,
                                               const Configuration &conf,
                                               std::string *ErrorStr) {
#ifdef LLVM_MODULE_MATERIALIZE_ALL_PERMANENTLY_ERRORCODE_BOOL
  if(std::error_code EC = M->materializeAllPermanently()){
    // We got an error, just return 0
    if(ErrorStr) *ErrorStr = EC.message();
    return 0;
  }
#elif defined LLVM_MODULE_MATERIALIZE_ALL_PERMANENTLY_BOOL_STRPTR
  if (M->MaterializeAllPermanently(ErrorStr)){
    // We got an error, just return 0
    return 0;
  }
#elif defined LLVM_MODULE_MATERIALIZE_LLVM_ALL_ERROR
  if (llvm::Error Err = M->materializeAll()) {
    std::string Msg;
    handleAllErrors(std::move(Err), [&](llvm::ErrorInfoBase &EIB) {
      Msg = EIB.message();
    });
    if (ErrorStr)
      *ErrorStr = Msg;
    // We got an error, just return 0
    return nullptr;
  }
#else
  if(std::error_code EC = M->materializeAll()){
    // We got an error, just return 0
    if(ErrorStr) *ErrorStr = EC.message();
    return 0;
  }
#endif

  return new ZInterpreterPSO(M,TB,conf);
}


void ZInterpreterPSO::runAux(int proc, int aux){
  /* Perform an update from store buffer to memory. */

  assert(0 <= aux && aux < int(pso_threads[proc].aux_to_byte.size()));

  SymAddr b0 = pso_threads[proc].aux_to_byte[aux];

  assert(pso_threads[proc].store_buffers.count(b0));
  assert(pso_threads[proc].store_buffers[b0].size());

  SymAddrSize ml = pso_threads[proc].store_buffers[b0].front().ml;
  assert(ml.addr == b0);

  // We consider all Stack writes as invisible and atomic,
  // hence we do not put them into the store buffer
  assert((b0.block.is_global() || b0.block.is_heap()) &&
         "No Stack writes are allowed in the store buffer");

  SymData sd(ml, ml.size);
  uint8_t *blk = (uint8_t*)sd.get_block();

  for(SymAddr b : ml){
    assert(pso_threads[proc].store_buffers.count(b));
    std::vector<PendingStoreByte> &sb = pso_threads[proc].store_buffers[b];
    assert(sb.size());
    assert(sb.front().ml == ml);
    blk[unsigned(b-b0)] = sb.front().val;
    assert(!DryRun); /**/
    if(!DryRun) {
      sb.erase(sb.begin());
      if(sb.empty()) pso_threads[proc].store_buffers.erase(b);
    }
  }

  TB.atomic_store(sd);
  assert(!DryRun); /**/

  if(!CheckedMemCpy(pso_threads[proc].aux_to_addr[aux],blk,ml.size)){
    llvm::errs() << "Interpreter: CheckedMemCpy failed during store-update\n";
    abort(); /**/
  }

  /* Should we reenable the thread after awaiting buffer flush? */
  switch(pso_threads[proc].awaiting_buffer_flush){
  case PSOThread::BFL_NO: // Do nothing
    break;
  case PSOThread::BFL_FULL:
    if(pso_threads[proc].all_buffers_empty()){
      /* Enable */
      pso_threads[proc].awaiting_buffer_flush = PSOThread::BFL_NO;
      TB.mark_available(proc);
    }
    break;
  case PSOThread::BFL_PARTIAL:
    if(ml.overlaps(pso_threads[proc].buffer_flush_ml) &&
       pso_threads[proc].readable(pso_threads[proc].buffer_flush_ml)){
      /* Enable */
      pso_threads[proc].awaiting_buffer_flush = PSOThread::BFL_NO;
      TB.mark_available(proc);
    }
    break;
  }

  if(pso_threads[proc].all_buffers_empty()){
    /* If the real thread has terminated, then wake up other threads
     * which are waiting to join with this one. */
    if(Threads[proc].ECStack.empty()){
      for(int p : Threads[proc].AwaitingJoin){
        TB.mark_available(p);
      }
    }
  }
}


bool ZInterpreterPSO::checkRefuse(llvm::Instruction &I) {
  int tid;
  if(isPthreadJoin(I,&tid)){
    if(0 <= tid && tid < int(Threads.size()) && tid != CurrentThread){
      if(Threads[tid].ECStack.size() ||
         !pso_threads[tid].all_buffers_empty() ||
         Threads[tid].AssumeBlocked){
        /* The awaited thread is still executing or assume-blocked. */
        TB.refuse_schedule();
        Threads[tid].AwaitingJoin.push_back(CurrentThread);
        return true;
      }
    }else{
      // Erroneous thread id
      // Allow execution (will produce an error trace)
    }
  }
  /* Refuse if I has fence semantics and the store buffer is
   * non-empty.
   */
  if(isFence(I) && !pso_threads[CurrentThread].all_buffers_empty()){
    pso_threads[CurrentThread].awaiting_buffer_flush = PSOThread::BFL_FULL;
    TB.refuse_schedule();
    return true;
  }
  /* Refuse if I is a load and the latest entry in the store buffer
   * which overlaps with the memory location targeted by I does not
   * target precisely the same bytes as I.
   */
  if(llvm::isa<llvm::LoadInst>(I)){
    llvm::ExecutionContext &SF = ECStack()->back();
    llvm::GenericValue SRC = getOperandValue(static_cast<llvm::LoadInst&>(I).getPointerOperand(), SF);
    llvm::GenericValue *Ptr = (llvm::GenericValue*)GVTOP(SRC);
    SymAddrSize mr = GetSymAddrSize(Ptr,static_cast<llvm::LoadInst&>(I).getType());
    if(!pso_threads[CurrentThread].readable(mr)){
      llvm::errs() << "Interpreter: Partial but not full overlap of a memory location for a store and load\n";
      abort(); /**/
      /* Block until this store buffer entry has disappeared from
       * the buffer.
       */
      pso_threads[CurrentThread].awaiting_buffer_flush = PSOThread::BFL_PARTIAL;
      pso_threads[CurrentThread].buffer_flush_ml = mr;
      TB.refuse_schedule();
      return true;
    }
  }
  /* Refuse if I is a lock and the mutex is already locked
   */
  llvm::GenericValue *ptr;
  if(isPthreadMutexLock(I,&ptr)){
    if(PthreadMutexes.count(ptr) &&
       PthreadMutexes[ptr].isLocked()){
      // Trying to access a locked mutex.
      TB.mutex_lock_fail({GetSymAddr(ptr),1});
      TB.refuse_schedule();
      //TB.mark_unavailable(CurrentThread);
      PthreadMutexes[ptr].waiting.insert(CurrentThread);
      return true;
    }else{
      // Either unlocked mutex, or uninitialized mutex.
      // In both cases let callPthreadMutex handle it.
    }
  }
  return false;
}


void ZInterpreterPSO::visitLoadInst(llvm::LoadInst &I){
  llvm::ExecutionContext &SF = ECStack()->back();
  llvm::GenericValue SRC = getOperandValue(I.getPointerOperand(), SF);
  llvm::GenericValue *Ptr = (llvm::GenericValue*)GVTOP(SRC);
  llvm::GenericValue Result;

  SymAddrSize ml = GetSymAddrSize(Ptr,I.getType());
  /**/ //TB.load(ml);

  assert(!DryRun); /**/

  /* Check store buffer for ROWE opportunity. */
  if(pso_threads[CurrentThread].store_buffers.count(ml.addr)){
    // We consider all Stack writes as invisible and atomic,
    // hence we do not put them into the store buffer
    assert(ml.addr.block.is_global() && ml.addr.block.is_heap());
    uint8_t *blk = new uint8_t[ml.size];
    for(SymAddr b : ml){
      assert(pso_threads[CurrentThread].store_buffers[b].back().ml == ml);
      blk[b - ml.addr] = pso_threads[CurrentThread].store_buffers[b].back().val;
    }
    CheckedLoadValueFromMemory(Result,(llvm::GenericValue*)blk,I.getType());
    SetValue(&I, Result, SF);
    // Loading value Result.IntVal.getSExtValue()
    TB.load(ml, (int) Result.IntVal.getSExtValue()); /**/
    delete[] blk;
    return;
  }

  /* Load from memory */
  bool res = CheckedLoadValueFromMemory(Result, Ptr, I.getType());
  if (!res) {
    llvm::errs() << "Interpreter: CheckedLoadValueFromMemory failed during load from memory\n";
    abort(); /**/
  }
  SetValue(&I, Result, SF);
  // Loading value Result.IntVal.getSExtValue()
  TB.load(ml, (int) Result.IntVal.getSExtValue()); /**/
}


void ZInterpreterPSO::visitStoreInst(llvm::StoreInst &I){
  llvm::ExecutionContext &SF = ECStack()->back();
  llvm::GenericValue Val = getOperandValue(I.getOperand(0), SF);
  llvm::GenericValue *Ptr = (llvm::GenericValue *)GVTOP
    (getOperandValue(I.getPointerOperand(), SF));

  SymData mb = GetSymData(Ptr, I.getOperand(0)->getType(), Val);

  const SymAddrSize& ml = mb.get_ref();
  // Stores to local memory on stack may not conflict

  PSOThread &thr = pso_threads[CurrentThread];

  if(I.getOrdering() == LLVM_ATOMIC_ORDERING_SCOPE::SequentiallyConsistent ||
     0 <= AtomicFunctionCall ||
     (!ml.addr.block.is_global() && !ml.addr.block.is_heap())) {
    /* Atomic store */
    /* We consider all Stack writes as invisible and atomic */
    if (ml.addr.block.is_global() || ml.addr.block.is_heap()) {
      llvm::errs() << "Interpreter: No support for visible Atomic store\n";
      abort(); /**/
    }
    //assert(thr.all_buffers_empty());
    assert(!DryRun); /**/
    CheckedStoreValueToMemory(Val, Ptr, I.getOperand(0)->getType());
  } else {
    /* Store to buffer */
    // Storing value Val.IntVal.getSExtValue()
    assert(ml.addr.block.is_global() || ml.addr.block.is_heap());
    if(thr.byte_to_aux.count(ml.addr) == 0){
      thr.byte_to_aux[ml.addr] = int(thr.aux_to_byte.size());
      thr.aux_to_byte.push_back(ml.addr);
      thr.aux_to_addr.push_back((uint8_t*)Ptr);
    }

    TB.store(mb, (int) Val.IntVal.getSExtValue());
    assert(!DryRun); /**/
    for(SymAddr b : ml){
      unsigned i = b - ml.addr;
      thr.store_buffers[b].push_back(PendingStoreByte(ml,((uint8_t*)mb.get_block())[i]));
    }
  }
}


void ZInterpreterPSO::visitFenceInst(llvm::FenceInst &I){
  if(I.getOrdering() == LLVM_ATOMIC_ORDERING_SCOPE::SequentiallyConsistent){
    TB.fence();
  }
}


void ZInterpreterPSO::visitInlineAsm(llvm::CallSite &CS, const std::string &asmstr){
  if(asmstr == "mfence"){
    TB.fence();
  }else if(asmstr == ""){ // Do nothing
  }else{
    throw std::logic_error("Unsupported inline assembly: "+asmstr);
  }
}
