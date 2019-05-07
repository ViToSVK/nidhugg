#ifndef __Z_HELPERS_H__
#define __Z_HELPERS_H__

#include <llvm/IR/Instructions.h>
#include <llvm/IR/Function.h>

#include "ZEvent.h"
#include "ZBasis.h"

/* helper functions */
inline bool sameMl(const ZEvent *ev1, const ZEvent *ev2) {
  assert(ev1 && ev2);
  return ev1->ml == ev2->ml;
}
inline bool sameMl(const ZEvent& ev1, const ZEvent& ev2) {
  return ev1.ml == ev2.ml;
}
inline bool sameMl(const Node *nd1, const Node *nd2) {
  assert(nd1 && nd2);
  return sameMl(nd1->getEvent(), nd2->getEvent());
}

//

inline bool isLoad(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::LOAD;
}
inline bool isLoad(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::LOAD;
}
inline bool isLoad(const Node *nd) {
  assert(nd);
  return isLoad(nd->getEvent());
}

//

inline bool isStore(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::STORE;
}
inline bool isStore(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::STORE;
}
inline bool isStore(const Node *nd) {
  assert(nd);
  return isStore(nd->getEvent());
}

//

inline bool isSpawn(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::SPAWN;
}
inline bool isSpawn(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::SPAWN;
}
inline bool isSpawn(const Node *nd) {
  assert(nd);
  return isSpawn(nd->getEvent());
}

//

inline bool isJoin(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::JOIN;
}
inline bool isJoin(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::JOIN;
}
inline bool isJoin(const Node *nd) {
  assert(nd);
  return isJoin(nd->getEvent());
}

//

inline bool isMutexInit(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_INIT;
}
inline bool isMutexInit(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_INIT;
}
inline bool isMutexInit(const Node *nd) {
  assert(nd);
  return isMutexInit(nd->getEvent());
}

//

inline bool isMutexDestroy(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_DESTROY;
}
inline bool isMutexDestroy(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_DESTROY;
}
inline bool isMutexDestroy(const Node *nd) {
  assert(nd);
  return isMutexDestroy(nd->getEvent());
}

//

inline bool isLock(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_LOCK;
}
inline bool isLock(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_LOCK;
}
inline bool isLock(const Node *nd) {
  assert(nd);
  return isLock(nd->getEvent());
}

//

inline bool isUnlock(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_UNLOCK;
}
inline bool isUnlock(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_UNLOCK;
}
inline bool isUnlock(const Node *nd) {
  assert(nd);
  return isUnlock(nd->getEvent());
}

//

inline bool isRead(const ZEvent *ev) {
  assert(ev);
  return isLoad(ev); // || isLock(ev) - treat them separately
}
inline bool isRead(const ZEvent& ev) {
  return isLoad(ev); // || isLock(ev) - treat them separately
}
inline bool isRead(const Node *nd) {
  assert(nd);
  return isRead(nd->getEvent());
}

//

inline bool isWrite(const ZEvent *ev) {
  assert(ev);
  return isStore(ev); //  || isUnlock(ev) - treat them separately
}
inline bool isWrite(const ZEvent& ev) {
  return isStore(ev); //  || isUnlock(ev) - treat them separately
}
inline bool isWrite(const Node *nd) {
  assert(nd);
  return isWrite(nd->getEvent());
}

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
