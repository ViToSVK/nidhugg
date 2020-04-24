#ifndef __Z_HELPERS_H__
#define __Z_HELPERS_H__

#include <llvm/IR/Instructions.h>
#include <llvm/IR/Function.h>

#include "ZEvent.h"


/* helper functions */

inline bool is_initial(const ZEvent *ev) {
  assert(ev);
  assert((ev->kind != ZEvent::Kind::INITIAL && ev->event_id() >= 0) ||
         (ev->value() == 0 && ev->event_id() == -1));
  return ev->kind == ZEvent::Kind::INITIAL;
}
inline bool is_initial(const ZEvent& ev) {
  assert((ev.kind != ZEvent::Kind::INITIAL && ev.event_id() >= 0) ||
         (ev.value() == 0 && ev.event_id() == -1));
  return ev.kind == ZEvent::Kind::INITIAL;
}

//

inline bool same_ml(const ZEvent *ev1, const ZEvent *ev2) {
  assert(ev1 && ev2);
  assert(!is_initial(ev1) && !is_initial(ev2));
  return ev1->ml() == ev2->ml();
}
inline bool same_ml(const ZEvent& ev1, const ZEvent& ev2) {
  assert(!is_initial(ev1) && !is_initial(ev2));
  return ev1.ml() == ev2.ml();
}

//

inline bool is_read(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::READ;
}
inline bool is_read(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::READ;
}

//

inline bool is_write(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::WRITE;
}
inline bool is_write(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::WRITE;
}

//

inline bool is_spawn(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::SPAWN;
}
inline bool is_spawn(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::SPAWN;
}

//

inline bool is_join(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::JOIN;
}
inline bool is_join(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::JOIN;
}

//

inline bool is_mutex_init(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_INIT;
}
inline bool is_mutex_init(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_INIT;
}

//

inline bool is_mutex_destroy(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_DESTROY;
}
inline bool is_mutex_destroy(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_DESTROY;
}

//

inline bool is_lock(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_LOCK;
}
inline bool is_lock(const ZEvent& ev) {
  return ev.kind == ZEvent::Kind::M_LOCK;
}

//

inline bool is_unlock(const ZEvent *ev) {
  assert(ev);
  return ev->kind == ZEvent::Kind::M_UNLOCK;
}
inline bool is_unlock(const ZEvent& ev) {
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

inline bool is_lock(const llvm::Instruction& I) {
  return is_function_call(&I, "pthread_mutex_lock");
}

inline bool is_unlock(const llvm::Instruction& I) {
  return is_function_call(&I, "pthread_mutex_unlock");
}

inline bool is_read(const llvm::Instruction& I) {
  return llvm::isa<llvm::LoadInst>(&I) || is_unlock(I);
}

inline bool is_write(const llvm::Instruction& I) {
  return llvm::isa<llvm::StoreInst>(&I) || is_lock(I);
}

inline bool is_global_load(const llvm::Instruction *I) {
  // XXX: what about pointers!
  return llvm::isa<llvm::LoadInst>(I) &&
          !llvm::isa<llvm::AllocaInst>(I->getOperand(0)->stripInBoundsOffsets());
}

inline bool is_join(const llvm::Instruction *I) {
  if (!I)
    return false;
  return is_function_call(I, "pthread_join");
}

inline bool is_pthread_create(const llvm::Instruction *I) {
  if (!I)
    return false;
  return is_function_call(I, "pthread_create");
}

#endif // __Z_HELPERS_H__
