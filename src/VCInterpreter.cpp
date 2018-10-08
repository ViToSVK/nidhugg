/* Copyright (C) 2014-2016 Carl Leonardsson
 * Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2018 Viktor Toman
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

#include<iostream>

#include "VCInterpreter.h"
#include "VCTraceBuilder.h"
#include "VCEvent.h"


static void SetValue(llvm::Value *V,
                     llvm::GenericValue Val, llvm::ExecutionContext &SF)
{
  SF.Values[V] = Val;
}

VCInterpreter::VCInterpreter(llvm::Module *M, VCTraceBuilder &TB,
                             const Configuration &conf)
  : TSOInterpreter(M,TB,conf), TB(TB)
{
}

llvm::ExecutionEngine *VCInterpreter::create(llvm::Module *M, VCTraceBuilder &TB,
                                              const Configuration &conf,
                                              std::string *ErrorStr)
{
#ifdef LLVM_MODULE_MATERIALIZE_ALL_PERMANENTLY_ERRORCODE_BOOL
  if(std::error_code EC = M->materializeAllPermanently()){
    // We got an error, just return 0
    if(ErrorStr) *ErrorStr = EC.message();
    return nullptr;
  }
#elif defined LLVM_MODULE_MATERIALIZE_ALL_PERMANENTLY_BOOL_STRPTR
  if (M->MaterializeAllPermanently(ErrorStr)){
    // We got an error, just return 0
    return nullptr;
  }
#elif LLVM_VERSION_MAJOR == 4
  if(llvm::Error EC = M->materializeAll()){
    // We got an error, just return 0
    if(EC)
      return 0;
  }
#else
  if(std::error_code EC = M->materializeAll()){
    // We got an error, just return 0
    if(ErrorStr) *ErrorStr = EC.message();
    return nullptr;
  }
#endif

  return new VCInterpreter(M,TB,conf);
}

void VCInterpreter::visitLoadInst(llvm::LoadInst &I)
{
  using namespace llvm;

  ExecutionContext &SF = ECStack()->back();
  GenericValue SRC = getOperandValue(I.getPointerOperand(), SF);
  GenericValue *Ptr = (GenericValue*)GVTOP(SRC);
  GenericValue Result;

  SymAddrSize Ptr_sas = GetSymAddrSize(Ptr,I.getType());

  if(!CheckedLoadValueFromMemory(Result, Ptr, I.getType()))
    return;

  // Loading value Result.IntVal.getSExtValue()
  TB.load(Ptr_sas, (int) Result.IntVal.getSExtValue());

  SetValue(&I, Result, SF);
}

void VCInterpreter::visitStoreInst(llvm::StoreInst &I)
{
  using namespace llvm;

  ExecutionContext &SF = ECStack()->back();
  GenericValue Val = getOperandValue(I.getOperand(0), SF);
  GenericValue *Ptr = (GenericValue *)GVTOP(getOperandValue(I.getPointerOperand(), SF));

  // Storing value Val.IntVal.getSExtValue()
  SymData sd = GetSymData(Ptr, I.getOperand(0)->getType(), Val);
  TB.atomic_store(sd, (int) Val.IntVal.getSExtValue());

  assert(tso_threads[CurrentThread].store_buffer.empty());

  CheckedStoreValueToMemory(Val, Ptr, I.getOperand(0)->getType());
}

bool VCInterpreter::checkRefuse(llvm::Instruction &I)
{
  int tid;
  if(isPthreadJoin(I,&tid)){
    if(0 <= tid && tid < int(Threads.size()) && tid != CurrentThread){
      if(Threads[tid].ECStack.size()){
        // The awaited thread is still executing.
        TB.mark_unavailable(CurrentThread);
        TB.refuse_schedule();
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

  llvm::GenericValue *ptr;
  if(isPthreadMutexLock(I,&ptr)){
    if(PthreadMutexes.count(ptr) &&
       PthreadMutexes[ptr].isLocked()){
      assert(false && "why am I here?");
      TB.mark_unavailable(CurrentThread);
      // TB.refuse_schedule();
      PthreadMutexes[ptr].waiting.insert(CurrentThread);
      return true;
    }else{
      // Either unlocked mutex, or uninitialized mutex.
      // In both cases let callPthreadMutex handle it.
    }
  }
  return false;
}
