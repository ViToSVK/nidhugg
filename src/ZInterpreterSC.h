/* Copyright (C) 2014-2016 Carl Leonardsson
 * Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2020 Viktor Toman
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

#ifndef __Z_INTERPRETER_SC_H__
#define __Z_INTERPRETER_SC_H__

#include <config.h>

#include "TSOInterpreter.h"
#include "ZBuilderSC.h"


/* A ZInterpreterSC is an interpreter running under the SC
 * semantics. The execution should be guided by scheduling from a
 * ZBuilderSC.
 */
class ZInterpreterSC : public TSOInterpreter {

  ZBuilderSC& TB;

 public:
  explicit ZInterpreterSC(llvm::Module *M, ZBuilderSC &TB,
                          const Configuration &conf = Configuration::default_conf);
  virtual ~ZInterpreterSC();

  static llvm::ExecutionEngine *create(llvm::Module *M, ZBuilderSC &TB,
                                       const Configuration &conf = Configuration::default_conf,
                                       std::string *ErrorStr = 0);

  virtual void visitLoadInst(llvm::LoadInst &I);
  virtual void visitStoreInst(llvm::StoreInst &I);
  virtual void visitFenceInst(llvm::FenceInst &I);
  virtual void visitAtomicCmpXchgInst(llvm::AtomicCmpXchgInst &I) {
    llvm::errs() << "Interpreter: No support for AtomicCmpXchg\n";
    abort();
  }
  virtual void visitAtomicRMWInst(llvm::AtomicRMWInst &I) {
    llvm::errs() << "Interpreter: No support for AtomicRMW\n";
    abort();
  }
  virtual void visitInlineAsm(llvm::CallSite &CS, const std::string &asmstr);
 protected:
  virtual void runAux(int proc, int aux) {
    llvm::errs() << "Interpreter: No auxiliary threads in SC\n";
    abort();
  }
  virtual bool checkRefuse(llvm::Instruction &I);
  /*
  virtual int newThread(const CPid &cpid);
  virtual bool isFence(llvm::Instruction &I);
  virtual void terminate(llvm::Type *RetTy, llvm::GenericValue Result);
  */
};


#endif // __Z_INTERPRETER_SC_H__
