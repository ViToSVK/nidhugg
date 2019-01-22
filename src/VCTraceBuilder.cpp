/* Copyright (C) 2016-2017 Marek Chalupa
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

#include <iostream>

#include "Debug.h"
#include "VCTrace.h"
#include "VCTraceBuilder.h"
#include "VCExplorer.h"

bool VCTraceBuilder::reset()
{
  if (this->has_error()) {
    this->error_trace = this->get_trace();
    return true;
  }

  // Add lock event for every thread ending with a failed mutex lock attempt
  add_failed_lock_attempts();

  // Construct the explorer with:
  // the initial trace, this original TB, the star_root_index
  VCExplorer explorer = VCExplorer(std::move(prefix),
                                   !somethingToAnnotate.empty(),
                                   *this,
                                   this->star_root_index,
                                   this->previous_mutation_process_first,
                                   this->root_before_nonroots);

  // Call the main method
  bool error = explorer.explore();

  // Print the result statistics
  explorer.print_stats();

  return error;
}

/* *************************** */
/* SCHEDULING                  */
/* *************************** */

bool VCTraceBuilder::schedule(int *proc, int *aux, int *alt, bool *dryrun)
{
  // Not using these arguments
  *dryrun = false;
  *alt = 0;
  *aux = -1;
  this->dryrun = false;

  if(sch_replay) {
    assert(!sch_initial && !sch_extend);
    return schedule_replay_trace(proc);
  }

  assert(!sch_replay && (sch_initial || sch_extend));

  return schedule_arbitrarily(proc);
}

bool VCTraceBuilder::schedule_replay_trace(int *proc)
{
  assert(!replay_trace.empty());
  assert(prefix_idx < (int) replay_trace.size());

  // prefix_idx is the index into prefix. Since we are up until
  // now only replaying the event sequence of replay_trace,
  // prefix_idx is also the index into replay_trace pointing
  // to the event we are currently replaying

  if (prefix_idx == -1 || replay_trace[prefix_idx].size == prefix[prefix_idx].size) {
    ++prefix_idx;
    if (prefix_idx == (int) replay_trace.size()) {
      // We are done replaying replay_trace,
      // so now we want to get an extension
      sch_replay = false;
      sch_extend = true;
      // Decrease back prefix_idx, so that
      // schedule_arbitrarily can access the last
      // event in prefix - it will increase the prefix
      // again if needed
      --prefix_idx;
      return schedule_arbitrarily(proc);
    }
    // Next instruction is the beginning of a new event
    assert(replay_trace.size() > prefix.size());
    assert(prefix_idx == (int) prefix.size());
    unsigned p = replay_trace[prefix_idx].iid.get_pid();
    // Mark that thread p owns the new event
    threads[p].event_indices.push_back(prefix_idx);
    // Create the new event
    assert(replay_trace[prefix_idx].cpid ==
           threads[p].cpid && "Inconsistent scheduling");
    // if above fails, search for p' tied to r_t[p_i].cpid,
    // scan r_t down and swap p with p' in all found events
    assert(replay_trace[prefix_idx].instruction_order ==
           threads[p].executed_instructions + 1 && "Inconsistent scheduling");
    assert(replay_trace[prefix_idx].event_order ==
           threads[p].executed_events && "Inconsistent scheduling");

    prefix.emplace_back(IID<IPid>(IPid(p),threads[p].last_event_index()),
                        threads[p].cpid,
                        threads[p].executed_instructions + 1, // +1 so that first will be 1
                        threads[p].executed_events, // so that first will be 0
                        prefix.size());
    // Mark that thread executes a new event
    ++threads[p].executed_events;
    ++threads[p].executed_instructions;
    assert(prefix.back().id == prefix.size() - 1);
  } else {
    // Next instruction is a continuation of the current event
    assert(replay_trace[prefix_idx].size > prefix[prefix_idx].size);
    // Increase the size of the current event, it will be scheduled
    ++prefix[prefix_idx].size;
    assert(replay_trace[prefix_idx].size >=
           prefix[prefix_idx].size && "Inconsistent scheduling");
    unsigned p = replay_trace[prefix_idx].iid.get_pid();
    ++threads[p].executed_instructions;
  }

  assert((unsigned) prefix_idx < replay_trace.size());
  unsigned p = replay_trace[prefix_idx].iid.get_pid();
  assert(threads[p].available);

  // Here used to be:
  // ++threads[p].clock[p];
  // After refactoring, it should be:
  // threads[p].event_indices.push_back(prefix_idx);
  // But I don't know why here, I would suspect to
  // run this only when we have a new event
  // Mark that thread p executes a new instruction
  // ++threads[p].executed_instructions;

  bool ret = schedule_thread(proc, p);
  assert(ret && "Bug in scheduling: could not reproduce a given replay trace");

  return ret;
}

bool VCTraceBuilder::schedule_arbitrarily(int *proc)
{
  assert(!sch_replay && (sch_initial || sch_extend));

  assert(threads.size() % 2 == 0);
  if (somethingToAnnotate.size() < (threads.size() / 2)) {
    // We prefer scheduling threads that have not yet
    // seen a new unannotated event; this improves
    // chances that new unannotated reads observe in
    // this trace a write that we will subsequently
    // include in our partial order
    for(unsigned p = 0; p < threads.size(); p += 2) {
      if (!somethingToAnnotate.count(p))
        if (schedule_thread(proc, p))
          return true;
    }
  }

  for(unsigned p = 0; p < threads.size(); p += 2) {
    if (schedule_thread(proc, p))
      return true;
  }

  // We did not schedule anything
  prefix.shrink_to_fit();
  return false;
}

// Schedule the next instruction to be the next instruction from thread p
bool VCTraceBuilder::schedule_thread(int *proc, unsigned p)
{
  if (threads[p].available && !threads[p].sleeping &&
      (conf.max_search_depth < 0 || threads[p].last_event_index() < conf.max_search_depth)) {

    if (sch_initial || sch_extend) {
      update_prefix(p);
    } else {
      assert(sch_replay);
    }

    // set the scheduled thread
    *proc = p/2;

    return true;
  }

  return false;
}

// Used only when we are not replaying replay_trace
void VCTraceBuilder::update_prefix(unsigned p)
{
  // Here used to be:
  // ++threads[p].clock[p];
  // After refactoring, it should be:
  // threads[p].event_indices.push_back(prefix_idx);
  // But I don't know why here, I would suspect to
  // run this only when we have a new event
  ++threads[p].executed_instructions;

  if (prefix_idx != -1 && (int) p == curnode().iid.get_pid() &&
      !curnode().may_conflict) {
    // 1) This will not be the very first instruction
    // 2) Context switch is not going to happen now
    // 3) In this event we have only invisible instructions so far
    // Because of the above, we extend the current event
    assert(prefix_idx == (int) (prefix.size() - 1));
    ++prefix[prefix_idx].size;
  } else {
    // Because one of 1)2)3) doesn't hold, we create a new event
    ++prefix_idx;
    // Mark that thread p owns the new event
    threads[p].event_indices.push_back(prefix_idx);
    // Create the new event
    prefix.emplace_back(IID<IPid>(IPid(p),threads[p].last_event_index()),
                        threads[p].cpid,
                        threads[p].executed_instructions, // first will be 1
                        threads[p].executed_events, // first will be 0
                        prefix.size());
    ++threads[p].executed_events;
    assert(prefix.back().id == prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }

  if (annotationDependantThread.count(p)) {
    // This node can not become part of our partial order
    // as the thread is created after a so-far unannotated node
    curnode().include_in_po = false;
  }
}

/* ******************************************************* */
/* ADDRESSING FEEDBACK FROM INTERPRETER AFTER SCHEDULING   */
/* ******************************************************* */

void VCTraceBuilder::refuse_schedule()
{
  // llvm::errs() << " REFUSESCH";
  assert(prefix_idx == int(prefix.size())-1);
  assert(!prefix.back().may_conflict);

  unsigned p = curnode().iid.get_pid();

  // Here used to be:
  // --threads[p].event_indices[p];

  --threads[p].executed_instructions; //

  if (curnode().size == 1) {
    // Refused instruction wanted to create a new event,
    // therefore this entire event needs to be deleted
    // Unmark ownership of this new event
    assert((int) threads[p].event_indices.back() == prefix_idx);
    threads[p].event_indices.pop_back();
    --threads[p].executed_events; //
    // Delete this new event
    prefix.pop_back();
    --prefix_idx;
  } else {
    // Refused instruction wanted to extend the current event
    --curnode().size;
  }

  assert(prefix_idx == int(prefix.size())-1);

  // Mark this thread as unavailable since its next instruction
  // is the one we tried to schedule now, and it got refused
  mark_unavailable(p/2, -1);
}

void VCTraceBuilder::mayConflict(const SymAddrSize *ml)
{
  auto& curn = curnode();
  curn.may_conflict = true;
  if (ml) curn.ml = *ml;

  #ifndef NDEBUG
  bool consistent = (!sch_replay ||
                     replay_trace[prefix_idx].kind == prefix[prefix_idx].kind);
  assert(consistent);
  #endif
}

void VCTraceBuilder::metadata(const llvm::MDNode *md)
{
  if (md) curnode().md = md;
}

IID<CPid> VCTraceBuilder::get_iid() const
{
  IPid pid = curnode().iid.get_pid();
  int idx = curnode().iid.get_index();
  return IID<CPid>(threads[pid].cpid,idx);
}

void VCTraceBuilder::spawn()
{
  //llvm::errs() << " SPAWN";
  assert(curnode().kind == VCEvent::Kind::DUMMY);
  curnode().kind = VCEvent::Kind::SPAWN;
  mayConflict();
  // store the CPid of the new thread
  IPid parent_ipid = curnode().iid.get_pid();
  CPid child_cpid = CPS.spawn(threads[parent_ipid].cpid);
  curnode().childs_cpid = child_cpid;

  threads.emplace_back(child_cpid, prefix_idx); // second arg was threads[parent_ipid].event_indices
  unsigned child_ipid = threads.size() - 1;
  threads.emplace_back(CPS.new_aux(child_cpid), prefix_idx); // second arg was threads[parent_ipid].event_indices
  threads.back().available = false; // Empty store buffer

  if (somethingToAnnotate.count(curnode().iid.get_pid()) ||
      annotationDependantThread.count(curnode().iid.get_pid())) {
    // The whole spawned thread can not become part of our partial order
    // as its creation depends on a so-far unannotated node
    annotationDependantThread.insert(child_ipid);
  }
}

void VCTraceBuilder::join(int tgt_proc)
{
  //llvm::errs() << " JOIN";
  assert(curnode().kind == VCEvent::Kind::DUMMY);
  // We make sure that every join is an event with exactly
  // one instruction - only the join instruction itself
  assert(prefix_idx == (int) (prefix.size() - 1));
  if (prefix[prefix_idx].size > 1) {
    // Some invisible instructions happened before the
    // join instruction, we split this into two events
    unsigned p = curnode().iid.get_pid();
    --prefix[prefix_idx].size; // Substract the join instruction
    // Create a new event with only the join instruction
    ++prefix_idx;
    // Mark that thread p owns the new event
    threads[p].event_indices.push_back(prefix_idx);
    // Create the new event
    prefix.emplace_back(IID<IPid>(IPid(p),threads[p].last_event_index()),
                        threads[p].cpid,
                        threads[p].executed_instructions, // first will be 1
                        threads[p].executed_events, // first will be 0
                        prefix.size());
    ++threads[p].executed_events;
    assert(prefix.back().id == prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }

  assert(curnode().kind == VCEvent::Kind::DUMMY);
  curnode().kind = VCEvent::Kind::JOIN;
  mayConflict();
  curnode().childs_cpid = threads[2*tgt_proc].cpid;

  if (somethingToAnnotate.count(2*tgt_proc) ||
      annotationDependantJoin.count(2*tgt_proc)) {
    // This join can not become part of our partial order
    // as its success depends on a so-far unannotated node
    curnode().include_in_po = false;
    annotationDependantJoin.insert(curnode().iid.get_pid());
  }
}

void VCTraceBuilder::atomic_store(const SymData &sd, int val)
{
  // Stores to local memory on stack may not conflict
  // Global stores happening with just one thread existing may not conflict
  // but we have to record them as such anyway to keep track of them
  const SymAddrSize &ml = sd.get_ref();
  if (ml.addr.block.is_global() || ml.addr.block.is_heap()) {
    //llvm::errs() << " STORE_" << ml.to_string() << " ";
    assert(curnode().kind == VCEvent::Kind::DUMMY);
    curnode().kind = VCEvent::Kind::STORE;
    mayConflict(&ml);
    curnode().value = val;
  }
}

void VCTraceBuilder::load(const SymAddrSize &ml, int val)
{
  // Loads from local memory on stack may not conflict
  // Also global loads happening with just one thread existing may not conflict
  if ((ml.addr.block.is_global() || ml.addr.block.is_heap()) &&
      threads.size() > 2) { // 2 because two entries for each thread
    //llvm::errs() << " LOAD_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << " ";
    assert(curnode().kind == VCEvent::Kind::DUMMY);
    curnode().kind = VCEvent::Kind::LOAD;
    mayConflict(&ml);
    curnode().value = val;
    if (!sch_replay)
      somethingToAnnotate.insert(curnode().iid.get_pid());
  }
}

void VCTraceBuilder::fence()
{
  assert(curnode().iid.get_pid() % 2 == 0);
}

void VCTraceBuilder::mutex_init(const SymAddrSize &ml)
{
  //llvm::errs() << " M_INIT_" << ml.to_string() << " ";
  fence();
  assert(!mutexes.count(ml.addr));
  assert(curnode().kind == VCEvent::Kind::DUMMY);
  curnode().kind = VCEvent::Kind::M_INIT;
  mayConflict(&ml);
  mutexes[ml.addr] = Mutex(prefix_idx);
  mutexes[ml.addr].value = 31337; // default-unlocked mutex
}

void VCTraceBuilder::mutex_destroy(const SymAddrSize &ml)
{
  //llvm::errs() << " M_DESTROY_" << ml.to_string() << " ";
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)) {
    // Assume static initialization
    mutexes[ml.addr] = Mutex();
  }
  assert(mutexes.count(ml.addr));
  assert(curnode().kind == VCEvent::Kind::DUMMY);
  curnode().kind = VCEvent::Kind::M_DESTROY;
  mayConflict(&ml);
  mutexes.erase(ml.addr);
}

void VCTraceBuilder::mutex_unlock(const SymAddrSize &ml)
{
  //llvm::errs() << " M_UNLOCK_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << " ";
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)) {
    // Assume static initialization
    mutexes[ml.addr] = Mutex();
  }
  assert(mutexes.count(ml.addr));

  Mutex &mutex = mutexes[ml.addr];
  assert(0 <= mutex.last_access);
  assert(mutex.locked && "Unlocked resource got unlocked again");
  assert(mutex.value == - 1 - (curnode().iid.get_pid() * 1000000)
         && "Unlocked by different process than the one that locked");

  assert(curnode().kind == VCEvent::Kind::DUMMY);
  curnode().kind = VCEvent::Kind::M_UNLOCK;
  mayConflict(&ml);
  curnode().value = ( curnode().iid.get_pid() * 1000000 )
                    + curnode().event_order; // WRITE
                    // mutex unlocked by this event (>=0)

  mutex.last_access = prefix_idx;
  mutex.locked = false;
  mutex.value = curnode().value; // WRITE
}

void VCTraceBuilder::mutex_lock(const SymAddrSize &ml)
{
  //llvm::errs() << " M_LOCK_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << "_";
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)) {
    // Assume static initialization
    mutexes[ml.addr] = Mutex();
  }
  assert(mutexes.count(ml.addr));
  Mutex &mutex = mutexes[ml.addr];

  assert(curnode().kind == VCEvent::Kind::DUMMY);
  // We make sure that every lock is an event with exactly
  // one instruction - only the lock instruction itself
  assert(prefix_idx == (int) (prefix.size() - 1));
  if (prefix[prefix_idx].size > 1) {
    // Some invisible instructions happened before the
    // lock instruction, we split this into two events
    unsigned p = curnode().iid.get_pid();
    --prefix[prefix_idx].size; // Substract the lock instruction
    // Create a new event with only the lock instruction
    ++prefix_idx;
    // Mark that thread p owns the new event
    threads[p].event_indices.push_back(prefix_idx);
    // Create the new event
    prefix.emplace_back(IID<IPid>(IPid(p),threads[p].last_event_index()),
                        threads[p].cpid,
                        threads[p].executed_instructions, // first will be 1
                        threads[p].executed_events, // first will be 0
                        prefix.size());
    ++threads[p].executed_events;
    assert(prefix.back().id == prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }

  assert(curnode().size == 1);
  assert(curnode().kind == VCEvent::Kind::DUMMY);
  curnode().kind = VCEvent::Kind::M_LOCK;

  mayConflict(&ml);
  curnode().value = mutex.value; // READ
  if (!sch_replay)
    somethingToAnnotate.insert(curnode().iid.get_pid());
  endsWithLockFail.erase(curnode().iid.get_pid());

  assert(!mutex.locked && mutex.value >= 0);
  mutex.last_lock = mutex.last_access = prefix_idx;
  mutex.locked = true;
  mutex.value = - 1 - (curnode().iid.get_pid() * 1000000);
                // mutex locked by this process (<0)

}

void VCTraceBuilder::mutex_lock_fail(const SymAddrSize &ml) {
  //llvm::errs() << " M_LOCKFAIL" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << "_";
  assert(!dryrun);
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)){
    // Assume static initialization
    mutexes[ml.addr] = Mutex();
  }
  assert(mutexes.count(ml.addr));
  #ifndef NDEBUG
  Mutex &mutex = mutexes[ml.addr];
  assert(0 <= mutex.last_lock);
  assert(mutex.locked && mutex.value < 0);
  #endif

  if (!sch_replay)
    somethingToAnnotate.insert(curnode().iid.get_pid());
  endsWithLockFail.emplace(curnode().iid.get_pid(), ml);
}

std::pair<std::vector<VCEvent>, bool> VCTraceBuilder::extendGivenTrace() {
  assert(sch_replay && !replay_trace.empty());

  std::unique_ptr<llvm::ExecutionEngine> EE(DPORDriver::create_execution_engine(M, *this, config));

  // Run main.
  EE->runFunctionAsMain(M->getFunction("main"), {"prog"}, 0);

  // Run static destructors.
  EE->runStaticConstructorsDestructors(true);

  // Add lock event for every thread ending with a failed mutex lock attempt
  add_failed_lock_attempts();

  return {prefix, !somethingToAnnotate.empty()};
}

void VCTraceBuilder::add_failed_lock_attempts() {
  // Add lock event for every thread ending with a failed mutex lock attempt
  for (auto p_ml : endsWithLockFail) {
    unsigned p = p_ml.first;
    ++threads[p].executed_instructions;
    ++prefix_idx;
    // Mark that thread p owns the new event
    threads[p].event_indices.push_back(prefix_idx);
    // Create the new event
    prefix.emplace_back(IID<IPid>(IPid(p),threads[p].last_event_index()),
                        threads[p].cpid,
                        threads[p].executed_instructions, // first will be 1
                        threads[p].executed_events, // first will be 0
                        prefix.size());
    ++threads[p].executed_events;
    assert(prefix.back().id == prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
    assert(curnode().size == 1);
    assert(curnode().kind == VCEvent::Kind::DUMMY);
    curnode().kind = VCEvent::Kind::M_LOCK;

    mayConflict(&(p_ml.second));
  }
}

Trace *VCTraceBuilder::get_trace() const
{
  if (error_trace) {
    assert(errors.size() == 0);
    return error_trace;
  }

  std::vector<IID<CPid> > cmp;
  std::vector<const llvm::MDNode*> cmp_md;
  std::vector<Error*> errs;
  if (errors.size() == 0)
    return nullptr;

  for(unsigned i = 0; i < prefix.size(); ++i) {
    cmp.push_back(IID<CPid>(threads[prefix[i].iid.get_pid()].cpid,
                            prefix[i].iid.get_index()));
    cmp_md.push_back(prefix[i].md);
  }
  for(unsigned i = 0; i < errors.size(); ++i){
    errs.push_back(errors[i]->clone());
  }
  Trace *t = new IIDSeqTrace(cmp,cmp_md,errs);
  t->set_blocked(false);
  return t;
}
