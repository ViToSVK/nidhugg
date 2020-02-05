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

#include "ZInterpreterSC.h"
#include "ZBuilderSC.h"


static void SetValue(llvm::Value *V, llvm::GenericValue Val, llvm::ExecutionContext &SF) {
  SF.Values[V] = Val;
}


ZInterpreterSC::ZInterpreterSC(llvm::Module *M, ZBuilderSC &TB,
                               const Configuration &conf)
  : TSOInterpreter(M,TB,conf), TB(TB) {}


ZInterpreterSC::~ZInterpreterSC() {}


llvm::ExecutionEngine *ZInterpreterSC::create(llvm::Module *M, ZBuilderSC &TB,
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

  return new ZInterpreterSC(M,TB,conf);
}


void ZInterpreterSC::runAux(int proc, int aux) {
  /* Perform an update from store buffer to memory. */
  assert(false && "Auxiliary threads should not exist (and be called) in SC");

  void *ref = tso_threads[proc].store_buffer.front().first;
  const SymData &blk = tso_threads[proc].store_buffer.front().second;
  TB.atomic_store(blk);

}


bool ZInterpreterSC::checkRefuse(llvm::Instruction &I) {
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
  /* Store buffer should always be empty in SC,
   * thus no need to check for fence-refuse.
   */
  assert(tso_threads[CurrentThread].store_buffer.empty() && "Empty buffer in SC");
  /* Partial but not full overlap of a memory location:
   * No need to check since buffer is always empty in SC.
   */
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


void ZInterpreterSC::visitLoadInst(llvm::LoadInst &I){
  llvm::ExecutionContext &SF = ECStack()->back();
  llvm::GenericValue SRC = getOperandValue(I.getPointerOperand(), SF);
  llvm::GenericValue *Ptr = (llvm::GenericValue*)GVTOP(SRC);
  llvm::GenericValue Result;

  SymAddrSize Ptr_sas = GetSymAddrSize(Ptr,I.getType());

  assert(!DryRun); /**/
  assert(tso_threads[CurrentThread].store_buffer.empty() && "Empty buffer in SC");

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


void ZInterpreterSC::visitStoreInst(llvm::StoreInst &I){
  llvm::ExecutionContext &SF = ECStack()->back();
  llvm::GenericValue Val = getOperandValue(I.getOperand(0), SF);
  llvm::GenericValue *Ptr = (llvm::GenericValue *)GVTOP
    (getOperandValue(I.getPointerOperand(), SF));

  Option<SymAddrSize> Ptr_sas = TryGetSymAddrSize(Ptr,I.getOperand(0)->getType());
  if (!Ptr_sas) {
    TB.segmentation_fault_error();
    abort();
    return;
  }

  assert(tso_threads[CurrentThread].store_buffer.empty() && "Empty buffer in SC");
  SymData sd = GetSymData(Ptr, I.getOperand(0)->getType(), Val);

  const SymAddrSize& ml = sd.get_ref();
  // Stores to local memory on stack may not conflict
  if(I.getOrdering() == LLVM_ATOMIC_ORDERING_SCOPE::SequentiallyConsistent ||
     0 <= AtomicFunctionCall ||
     (!ml.addr.block.is_global() && !ml.addr.block.is_heap())) {
    /* Stack write */
    if (ml.addr.block.is_global() || ml.addr.block.is_heap()) {
      llvm::errs() << "Interpreter: No support for visible Atomic store\n";
      abort(); /**/
    }
    //assert(tso_threads[CurrentThread].store_buffer.empty());
    assert(!DryRun); /**/
    CheckedStoreValueToMemory(Val, Ptr, I.getOperand(0)->getType());
  } else {
    /* Heap or Global write */
    // Storing value Val.IntVal.getSExtValue()
    assert(ml.addr.block.is_global() || ml.addr.block.is_heap());
    TB.store(sd, (int) Val.IntVal.getSExtValue()); // Calls TB.atomic_store

    assert(!DryRun); /**/
    CheckedStoreValueToMemory(Val, Ptr, I.getOperand(0)->getType());
  }
}


void ZInterpreterSC::visitFenceInst(llvm::FenceInst &I){
  if(I.getOrdering() == LLVM_ATOMIC_ORDERING_SCOPE::SequentiallyConsistent){
    TB.fence();
  }
}


void ZInterpreterSC::visitInlineAsm(llvm::CallSite &CS, const std::string &asmstr){
  if(asmstr == "mfence"){
    TB.fence();
  }else if(asmstr == ""){ // Do nothing
  }else{
    throw std::logic_error("Unsupported inline assembly: "+asmstr);
  }
}
