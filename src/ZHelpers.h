#ifndef __Z_HELPERS_H__
#define __Z_HELPERS_H__

#include <llvm/IR/Instructions.h>
#include <llvm/IR/Function.h>

#include "ZEvent.h"

/* helper functions */
inline bool sameMl(const ZEvent *ev1, const ZEvent *ev2) {
  assert(ev1 && ev2);
  return ev1->ml == ev2->ml;
}
inline bool sameMl(const ZEvent& ev1, const ZEvent& ev2) {
  return ev1.ml == ev2.ml;
}

//

inline bool isRead(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::READ;
}
inline bool isRead(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::READ;
}

//

inline bool isWriteBuf(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::WRITEB;
}
inline bool isWriteBuf(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::WRITEB;
}


//

inline bool isWriteMem(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::WRITEM;
}
inline bool isWriteMem(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::WRITEM;
}


//

inline bool isSpawn(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::SPAWN;
}
inline bool isSpawn(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::SPAWN;
}

//

inline bool isJoin(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::JOIN;
}
inline bool isJoin(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::JOIN;
}

//

inline bool isMutexInit(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_INIT;
}
inline bool isMutexInit(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_INIT;
}

//

inline bool isMutexDestroy(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_DESTROY;
}
inline bool isMutexDestroy(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_DESTROY;
}

//

inline bool isLock(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_LOCK;
}
inline bool isLock(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_LOCK;
}

//

inline bool isUnlock(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_UNLOCK;
}
inline bool isUnlock(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_UNLOCK;
}

//

/* llvm instruction helpers */
inline bool is_function_call(const llvm::Instruction *I, const char *name) {
  using namespace llvm;

  const CallInst *CI = dyn_cast<CallInst>(I);
  if(!CI)
    return false;

  const Function *F
    = dyn_cast<Function>(CI->getCalledValue()->stripPointerCasts());
  if(!F || !F->getName().equals(name))
    return false;

  return true;
}

inline bool isLock(const llvm::Instruction& I) {
  return is_function_call(&I, "pthread_mutex_lock");
}

inline bool isUnlock(const llvm::Instruction& I) {
  return is_function_call(&I, "pthread_mutex_unlock");
}

inline bool isRead(const llvm::Instruction& I) {
  return llvm::isa<llvm::LoadInst>(&I) || isUnlock(I);
}

inline bool isWrite(const llvm::Instruction& I) {
  return llvm::isa<llvm::StoreInst>(&I) || isLock(I);
}

inline bool isGlobalLoad(const llvm::Instruction *I) {
  // XXX: what about pointers!
  return llvm::isa<llvm::LoadInst>(I) &&
          !llvm::isa<llvm::AllocaInst>(I->getOperand(0)->stripInBoundsOffsets());
}

inline bool isJoin(const llvm::Instruction *I) {
  if (!I)
    return false;
  return is_function_call(I, "pthread_join");
}

inline bool isPthreadCreate(const llvm::Instruction *I) {
  if (!I)
    return false;
  return is_function_call(I, "pthread_create");
}

#endif // __Z_HELPERS_H__
