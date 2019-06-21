/* Copyright (C) 2014-2016 Carl Leonardsson
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

#include "ZInterpreterTSO.h"
#include "ZBuilderTSO.h"


static void SetValue(llvm::Value *V, llvm::GenericValue Val, llvm::ExecutionContext &SF) {
  SF.Values[V] = Val;
}


ZInterpreterTSO::ZInterpreterTSO(llvm::Module *M, ZBuilderTSO &TB,
                                 const Configuration &conf)
  : TSOInterpreter(M,TB,conf), TB(TB) {}


ZInterpreterTSO::~ZInterpreterTSO() {}


llvm::ExecutionEngine *ZInterpreterTSO::create(llvm::Module *M, ZBuilderTSO &TB,
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

  return new ZInterpreterTSO(M,TB,conf);
}


void ZInterpreterTSO::runAux(int proc, int aux) {
  /* Perform an update from store buffer to memory. */

  assert(aux == 0);
  assert(tso_threads[proc].store_buffer.size());

  void *ref = tso_threads[proc].store_buffer.front().first;
  const SymData &blk = tso_threads[proc].store_buffer.front().second;

  // We consider all Stack writes as invisible and atomic,
  // hence we do not put them into the store buffer
  assert((blk.get_ref().addr.block.is_global() ||
          blk.get_ref().addr.block.is_heap()) &&
         "No Stack writes are allowed in the store buffer");

  TB.atomic_store(blk);

  assert(!DryRun); /**/

  if(!CheckedMemCpy((uint8_t*)ref,(uint8_t*)blk.get_block(),blk.get_ref().size)) {
    llvm::errs() << "Interpreter: CheckedMemCpy failed during store-update\n";
    abort(); /**/
  }

  for(unsigned i = 0; i < tso_threads[proc].store_buffer.size()-1; ++i){
    tso_threads[proc].store_buffer[i] = tso_threads[proc].store_buffer[i+1];
  }
  tso_threads[proc].store_buffer.pop_back();

  if(int(tso_threads[proc].store_buffer.size()) <= tso_threads[proc].partial_buffer_flush){
    assert(0 <= tso_threads[proc].partial_buffer_flush);
    /* The real thread was waiting for the buffer to flush. Enable the
     * real thread. */
    tso_threads[proc].partial_buffer_flush = -1;
    TB.mark_available(proc);
  }

  if(tso_threads[proc].store_buffer.empty()){
    /* Disable update thread. */
    TB.mark_unavailable(proc,0);
    /* If the real thread has terminated, then wake up other threads
     * which are waiting to join with this one. */
    if(Threads[proc].ECStack.empty()){
      for(int p : Threads[proc].AwaitingJoin){
        TB.mark_available(p);
      }
    }
  }
}


bool ZInterpreterTSO::checkRefuse(llvm::Instruction &I) {
  int tid;
  if(isPthreadJoin(I,&tid)){
    if(0 <= tid && tid < int(Threads.size()) && tid != CurrentThread){
      if(Threads[tid].ECStack.size() ||
         tso_threads[tid].store_buffer.size() ||
         Threads[tid].AssumeBlocked){
        /* The awaited thread is still executing or assume-blocked. */
        TB.refuse_schedule();
        // TB.mark_unavailable(CurrentThread);
        // Upon the threads termination, this thread is made
        // available in Interperter::terminate()
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
  if(isFence(I) && !tso_threads[CurrentThread].store_buffer.empty()){
    tso_threads[CurrentThread].partial_buffer_flush = 0;
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
    for(int i = int(tso_threads[CurrentThread].store_buffer.size())-1; 0 <= i; --i){
      if(mr.overlaps(tso_threads[CurrentThread].store_buffer[i].second.get_ref())){
        if(mr != tso_threads[CurrentThread].store_buffer[i].second.get_ref()){
          llvm::errs() << "Interpreter: Partial but not full overlap of a memory location for a store and load\n";
          abort(); /**/
          /* Block until this store buffer entry has disappeared from
           * the buffer.
           */
          tso_threads[CurrentThread].partial_buffer_flush =
            int(tso_threads[CurrentThread].store_buffer.size()) - i - 1;
          TB.refuse_schedule();
          return true;
        }
        break;
      }
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


void ZInterpreterTSO::visitLoadInst(llvm::LoadInst &I){
  llvm::ExecutionContext &SF = ECStack()->back();
  llvm::GenericValue SRC = getOperandValue(I.getPointerOperand(), SF);
  llvm::GenericValue *Ptr = (llvm::GenericValue*)GVTOP(SRC);
  llvm::GenericValue Result;

  SymAddrSize Ptr_sas = GetSymAddrSize(Ptr,I.getType());

  assert(!DryRun); /**/

  /* Check store buffer for ROWE opportunity. */
  for(int i = int(tso_threads[CurrentThread].store_buffer.size())-1; 0 <= i; --i){
    if(Ptr_sas.addr == tso_threads[CurrentThread].store_buffer[i].second.get_ref().addr){
      /* Read-Own-Write-Early */
      assert(GetSymAddrSize(Ptr,I.getType()).size == tso_threads[CurrentThread].store_buffer[i].second.get_ref().size);
      CheckedLoadValueFromMemory(Result,(llvm::GenericValue*)tso_threads[CurrentThread].store_buffer[i].second.get_block(),I.getType());
      SetValue(&I, Result, SF);
      // Loading value Result.IntVal.getSExtValue()
      TB.load(Ptr_sas, (int) Result.IntVal.getSExtValue()); /**/
      return;
    }
  }

  /* Load from memory */
  bool res = CheckedLoadValueFromMemory(Result, Ptr, I.getType());
  if (!res) {
    llvm::errs() << "Interpreter: CheckedLoadValueFromMemory failed during load from memory\n";
    abort(); /**/
  }
  SetValue(&I, Result, SF);
  // Loading value Result.IntVal.getSExtValue()
  TB.load(Ptr_sas, (int) Result.IntVal.getSExtValue()); /**/
}


void ZInterpreterTSO::visitStoreInst(llvm::StoreInst &I){
  llvm::ExecutionContext &SF = ECStack()->back();
  llvm::GenericValue Val = getOperandValue(I.getOperand(0), SF);
  llvm::GenericValue *Ptr = (llvm::GenericValue *)GVTOP
    (getOperandValue(I.getPointerOperand(), SF));

  SymData sd = GetSymData(Ptr, I.getOperand(0)->getType(), Val);

  const SymAddrSize& ml = sd.get_ref();
  // Stores to local memory on stack may not conflict
  if(I.getOrdering() == LLVM_ATOMIC_ORDERING_SCOPE::SequentiallyConsistent ||
     0 <= AtomicFunctionCall ||
     (!ml.addr.block.is_global() && !ml.addr.block.is_heap())) {
    /* Atomic store */
    /* We consider all Stack writes as invisible and atomic */
    if (ml.addr.block.is_global() || ml.addr.block.is_heap()) {
      llvm::errs() << "Interpreter: No support for visible Atomic store\n";
      abort(); /**/
    }
    //assert(tso_threads[CurrentThread].store_buffer.empty());
    assert(!DryRun); /**/
    CheckedStoreValueToMemory(Val, Ptr, I.getOperand(0)->getType());
  } else {
    /* Store to buffer */
    // Storing value Val.IntVal.getSExtValue()
    assert(ml.addr.block.is_global() || ml.addr.block.is_heap());
    TB.store(sd, (int) Val.IntVal.getSExtValue());
    assert(!DryRun); /**/
    tso_threads[CurrentThread].store_buffer.emplace_back(Ptr, std::move(sd));
  }
}


void ZInterpreterTSO::visitFenceInst(llvm::FenceInst &I){
  if(I.getOrdering() == LLVM_ATOMIC_ORDERING_SCOPE::SequentiallyConsistent){
    TB.fence();
  }
}


void ZInterpreterTSO::visitInlineAsm(llvm::CallSite &CS, const std::string &asmstr){
  if(asmstr == "mfence"){
    TB.fence();
  }else if(asmstr == ""){ // Do nothing
  }else{
    throw std::logic_error("Unsupported inline assembly: "+asmstr);
  }
}
