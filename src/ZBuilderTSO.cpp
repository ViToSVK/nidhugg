/* Copyright (C) 2016-2017 Marek Chalupa
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

#include "Debug.h"
#include "ZTrace.h"
#include "ZBuilderTSO.h"
#include "ZExplorer.h"


bool ZBuilderTSO::reset()
{
  if (this->has_error()) {
    this->error_trace = this->get_trace();
    return true;
  }

  // Add lock event for every thread ending with a failed mutex lock attempt
  add_failed_lock_attempts();

  // Construct the explorer with:
  // the initial trace, this original TB, the star_root_index
  ZExplorer explorer = ZExplorer(std::move(prefix),
                                 !somethingToAnnotate.empty(),
                                 *this,
                                 this->star_root_index);

  // Call the main method
  bool error = explorer.explore();

  // Print the result statistics
  // explorer.print_stats();

  return error;
}


/* *************************** */
/* SCHEDULING                  */
/* *************************** */


bool ZBuilderTSO::schedule(int *proc, int *aux, int *alt, bool *dryrun)
{
  // Not using these arguments
  *dryrun = false;
  *alt = 0;
  this->dryrun = false;

  if(sch_replay) {
    assert(!sch_extend);
    return schedule_replay_trace(proc, aux);
  }

  assert(!sch_replay && sch_extend);
  return schedule_arbitrarily(proc, aux);
}


bool ZBuilderTSO::schedule_replay_trace(int *proc, int *aux)
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
      return schedule_arbitrarily(proc, aux);
    }
    // Next instruction is the beginning of a new event
    assert(replay_trace.size() > prefix.size());
    assert(prefix_idx == (int) prefix.size());
    unsigned p = replay_trace[prefix_idx].iid.get_pid();
    // If below assertion fails, search for p' tied to r_t[p_i].cpid,
    // scan r_t down and swap p with p' in all found events
    assert(replay_trace[prefix_idx].cpid ==
           threads[p].cpid && "IPID<->CPID correspondence has changed");
    // Mark that thread p owns the new event
    threads[p].event_indices.push_back(prefix_idx);
    // Create the new event
    assert(replay_trace[prefix_idx].instruction_order ==
           threads[p].executed_instructions + 1 && "Inconsistent scheduling");
    assert(replay_trace[prefix_idx].eventID() ==
           threads[p].executed_events && "Inconsistent scheduling");

    prefix.emplace_back(IID<IPid>(IPid(p),threads[p].last_event_index()),
                        threads[p].cpid,
                        threads[p].executed_instructions + 1, // +1 so that first will be 1
                        threads[p].executed_events, // so that first will be 0
                        prefix.size());
    // Mark that thread executes a new event
    ++threads[p].executed_events;
    ++threads[p].executed_instructions;
    assert(prefix.back().traceID() == (int) prefix.size() - 1);
  } else {
    // Next instruction is a continuation of the current event
    assert(replay_trace[prefix_idx].size > prefix[prefix_idx].size);
    // Increase the size of the current event, it will be scheduled
    ++prefix[prefix_idx].size;
    assert(replay_trace[prefix_idx].size >=
           prefix[prefix_idx].size && "Inconsistent scheduling");
    unsigned p = replay_trace[prefix_idx].iid.get_pid();
    // Handle when next instruction is an invisible one
    // coming from an auxiliary thread
    if (p % 2 == 0 &&
        replay_trace[prefix_idx].aux_invisible.count(prefix[prefix_idx].size)) {
      p++;
      assert(p == replay_trace[prefix_idx].aux_invisible[prefix[prefix_idx].size]);
    }
    ++threads[p].executed_instructions;
  }

  assert((unsigned) prefix_idx < replay_trace.size());
  unsigned p = replay_trace[prefix_idx].iid.get_pid();
  // Handle when next instruction is an invisible one
  // coming from an auxiliary thread
  if (p % 2 == 0 &&
      replay_trace[prefix_idx].aux_invisible.count(prefix[prefix_idx].size)) {
    p++;
    assert(p == replay_trace[prefix_idx].aux_invisible[prefix[prefix_idx].size]);
  }
  assert(threads[p].available);

  // Here used to be:
  // ++threads[p].clock[p];
  // After refactoring, it should be:
  // threads[p].event_indices.push_back(prefix_idx);
  // But I don't know why here, I would suspect to
  // run this only when we have a new event
  // Mark that thread p executes a new instruction
  // ++threads[p].executed_instructions;

  bool ret = schedule_thread(proc, aux, p);
  assert(ret && "Bug in scheduling: could not reproduce a given replay trace");

  return ret;
}


bool ZBuilderTSO::schedule_arbitrarily(int *proc, int *aux)
{
  assert(!sch_replay && sch_extend);

  assert(threads.size() % 2 == 0);
  // We prefer scheduling threads that have not yet
  // seen a new unannotated event; this improves
  // chances that new unannotated reads observe in
  // this trace a write that we will subsequently
  // include in our partial order
  // Further we prefer auxiliary before real threads

  const unsigned sz = threads.size();
  unsigned p;
  if (somethingToAnnotate.size() < (threads.size() / 2)) {
    for (p = 1; p < sz; p += 2) { // Loop through auxiliary threads
      if (!somethingToAnnotate.count(p - 1))
        if (schedule_thread(proc, aux, p))
          return true;
    }
    for (p = 0; p < sz; p += 2) { // Loop through real threads
      if (!somethingToAnnotate.count(p))
        if (schedule_thread(proc, aux, p))
          return true;
    }
  }

  for (p = 1; p < sz; p += 2) { // Loop through auxiliary threads
    if (schedule_thread(proc, aux, p))
      return true;
  }

  for (p = 0; p < sz; p += 2) { // Loop through real threads
    if (schedule_thread(proc, aux, p))
      return true;
  }

  // We did not schedule anything
  prefix.shrink_to_fit();
  return false;
}


// Schedule the next instruction to be the next instruction from thread p
bool ZBuilderTSO::schedule_thread(int *proc, int *aux, unsigned p)
{
  if (threads[p].available && !threads[p].sleeping &&
      (conf.max_search_depth < 0 || threads[p].last_event_index() < conf.max_search_depth)) {

    if (sch_extend) {
      update_prefix(p);
    } else {
      assert(sch_replay);
    }

    // set the scheduled thread and aux
    *proc = p/2;
    *aux = (p % 2) - 1; // -1 for real, 0 for auxiliary

    return true;
  }

  return false;
}


// Used only when we are not replaying replay_trace
void ZBuilderTSO::update_prefix(unsigned p)
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
    assert(prefix.back().traceID() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }
}


/* ******************************************************* */
/* ADDRESSING FEEDBACK FROM INTERPRETER AFTER SCHEDULING   */
/* ******************************************************* */


void ZBuilderTSO::refuse_schedule()
{
  //llvm::errs() << " REFUSESCH";
  assert(prefix_idx == int(prefix.size())-1);
  assert(!prefix.back().may_conflict);

  unsigned p = curnode().iid.get_pid();
  assert(p % 2 == 0 && "Only real threads can be refused");

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


void ZBuilderTSO::mayConflict(const SymAddrSize *ml)
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


void ZBuilderTSO::metadata(const llvm::MDNode *md)
{
  if (md) curnode().md = md;
}


IID<CPid> ZBuilderTSO::get_iid() const
{
  IPid pid = curnode().iid.get_pid();
  int idx = curnode().iid.get_index();
  return IID<CPid>(threads[pid].cpid,idx);
}


void ZBuilderTSO::spawn()
{
  //llvm::errs() << " SPAWN";
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::SPAWN;
  mayConflict();
  // store the CPid of the new thread
  IPid parent_ipid = curnode().iid.get_pid();
  CPid child_cpid = CPS.spawn(threads[parent_ipid].cpid);
  curnode().childs_cpid = child_cpid;

  threads.emplace_back(child_cpid, prefix_idx); // second arg was threads[parent_ipid].event_indices
  threads.emplace_back(CPS.new_aux(child_cpid), prefix_idx); // second arg was threads[parent_ipid].event_indices
  threads.back().available = false; // New thread starts with an empty store buffer
}


void ZBuilderTSO::join(int tgt_proc)
{
  //llvm::errs() << " JOIN";
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
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
    assert(prefix.back().traceID() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }

  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::JOIN;
  mayConflict();
  curnode().childs_cpid = threads[2*tgt_proc].cpid;
}


void ZBuilderTSO::store(const SymData &sd, int val)
{
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);

  unsigned p = curnode().iid.get_pid();
  assert(p % 2 == 0);

  // Put the store to buffer
  threads[p].store_buffer.push_back(PendingStore(sd.get_ref(),prefix_idx,last_md));
  threads[p+1].available = true;

  const SymAddrSize& ml = sd.get_ref();
  // Stores to local memory on stack may not conflict
  // Global stores happening with just one thread existing may not conflict
  // but we have to record them as such anyway to keep track of them
  if (ml.addr.block.is_global() || ml.addr.block.is_heap()) {
    //llvm::errs() << " BUFFERWRITE_" << ml.to_string() << " ";
    assert(curnode().kind == ZEvent::Kind::DUMMY);
    curnode().kind = ZEvent::Kind::WRITEB;
    mayConflict(&ml);
    curnode().value = val;

    // Enqueue the update into visible store queue
    if (!visibleStoreQueue.count(p))
      visibleStoreQueue.emplace(p, std::vector<int>());
    visibleStoreQueue[p].push_back(prefix_idx);
  }
}


void ZBuilderTSO::atomic_store(const SymData &sd)
{
  assert(!dryrun);
  assert((curnode().iid.get_pid() % 2 == 1 ||
          (sch_replay && !sch_extend &&
           replay_trace[prefix_idx].aux_invisible.count(curnode().size))) &&
         "Only auxiliary threads can perform memory-writes");

  unsigned auxp = curnode().iid.get_pid();
  if (auxp % 2 != 1)
    auxp = replay_trace[prefix_idx].aux_invisible[curnode().size];
  assert(auxp % 2 == 1);
  unsigned realp = auxp - 1;
  assert(auxp > realp && auxp == realp + 1);

  // Remove pending store from buffer
  for(unsigned i=0; i<threads[realp].store_buffer.size()-1; ++i) {
    threads[realp].store_buffer[i] = threads[realp].store_buffer[i+1];
  }
  threads[realp].store_buffer.pop_back();
  if(threads[realp].store_buffer.empty()) {
    threads[auxp].available = false;
  }

  const SymAddrSize& ml = sd.get_ref();
  if (ml.addr.block.is_global() || ml.addr.block.is_heap()) {
    //llvm::errs() << " MEMORYWRITE_" << ml.to_string() << " ";
    assert(curnode().kind == ZEvent::Kind::DUMMY);
    curnode().kind = ZEvent::Kind::WRITEM;
    mayConflict(&ml);

    // Locate corresponding buffer-write
    assert(visibleStoreQueue.count(realp) && !visibleStoreQueue[realp].empty());
    curnode().write_other_trace_id = visibleStoreQueue[realp].front();
    assert(curnode().write_other_trace_id >= 0 &&
           curnode().write_other_trace_id <= (int) prefix.size() - 2);
    assert(prefix[curnode().write_other_trace_id].kind == ZEvent::Kind::WRITEB);
    assert(prefix[curnode().write_other_trace_id].ml == ml &&
           "Inconsistent memory locations");
    curnode().value = prefix[curnode().write_other_trace_id].value;
    assert(prefix[curnode().write_other_trace_id].write_other_trace_id == -1);
    prefix[curnode().write_other_trace_id].write_other_trace_id = prefix_idx;

    // Dequeue oldest update from queue
    for(unsigned i=0; i<visibleStoreQueue[realp].size()-1; ++i){
      visibleStoreQueue[realp][i] = visibleStoreQueue[realp][i+1];
    }
    visibleStoreQueue[realp].pop_back();

    lastWrite[ml] = prefix_idx;
  } else if (sch_extend) {
    assert(curnode().iid.get_pid() == (int) auxp);
    if (curnode().size == 1 && !curnode().may_conflict && prefix_idx > 0 &&
        prefix[prefix_idx - 1].iid.get_pid() == (int) realp &&
        !prefix[prefix_idx - 1].may_conflict &&
        prefix[prefix_idx - 1].kind == ZEvent::Kind::DUMMY) {
      // This is a single invisible auxiliary instruction,
      // previous event is invisible and belongs to
      // the corresponding real thread
      // We squash these two events together

      // Do not do this: --threads[auxp].executed_instructions;
      // Unmark ownership of this new event
      assert((int) threads[auxp].event_indices.back() == prefix_idx);
      threads[auxp].event_indices.pop_back();
      --threads[auxp].executed_events; //
      // Delete this new event
      prefix.pop_back();
      --prefix_idx;

      assert(curnode().iid.get_pid() == (int) realp &&
             !curnode().may_conflict &&
             curnode().kind == ZEvent::Kind::DUMMY);
      // We extend the now-current event
      assert(prefix_idx == (int) (prefix.size() - 1));
      ++prefix[prefix_idx].size;
      // We add that this instruction actually belongs
      // to its auxiliary thread, to keep track of it
      assert(!curnode().aux_invisible.count(curnode().size));
      curnode().aux_invisible.emplace(curnode().size, auxp);
    }
  }
}


void ZBuilderTSO::load(const SymAddrSize &ml, int val)
{
  // Loads from local memory on stack may not conflict
  // Also global loads happening with just one thread existing may not conflict
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  if ((ml.addr.block.is_global() || ml.addr.block.is_heap()) &&
      threads.size() > 2) { // 2 because two entries for each thread
    //llvm::errs() << " READ_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << " ";
    assert(curnode().kind == ZEvent::Kind::DUMMY);
    curnode().kind = ZEvent::Kind::READ;
    mayConflict(&ml);
    curnode().value = val;

    // Set ID of observed event
    // last memory-write
    int obs_idx = (lastWrite.count(ml)) ? lastWrite[ml] : -1;
    // check for buffer-writes present in queue
    unsigned p = curnode().iid.get_pid();
    if (!visibleStoreQueue.count(p))
      visibleStoreQueue.emplace(p, std::vector<int>());
    else {
      for(int i=visibleStoreQueue[p].size()-1; i>=0; --i) {
        int idx = visibleStoreQueue[p][i];
        if (prefix[idx].ml == ml) {
          obs_idx = idx;
          break;
        }
      }
    }
    curnode().observed_trace_id = obs_idx;

    if (!sch_replay)
      somethingToAnnotate.insert(curnode().iid.get_pid());
  }
}


void ZBuilderTSO::fence()
{
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  assert(!visibleStoreQueue.count(curnode().iid.get_pid()) ||
         visibleStoreQueue[curnode().iid.get_pid()].empty());
  assert(threads[curnode().iid.get_pid()].store_buffer.empty());
  curnode().fence = true;
}


void ZBuilderTSO::mutex_init(const SymAddrSize &ml)
{
  //llvm::errs() << " M_INIT_" << ml.to_string() << " ";
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  fence();
  assert(!mutexes.count(ml.addr));
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::M_INIT;
  mayConflict(&ml);
  mutexes[ml.addr] = Mutex(-1); // prefix_idx
}


void ZBuilderTSO::mutex_destroy(const SymAddrSize &ml)
{
  //llvm::errs() << " M_DESTROY_" << ml.to_string() << " ";
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)) {
    // Assume static initialization
    mutexes[ml.addr] = Mutex(-1);
  }
  assert(mutexes.count(ml.addr));
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::M_DESTROY;
  mayConflict(&ml);
  mutexes.erase(ml.addr);
}


void ZBuilderTSO::mutex_unlock(const SymAddrSize &ml)
{
  //llvm::errs() << " M_UNLOCK_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << " ";
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)) {
    // Assume static initialization
    mutexes[ml.addr] = Mutex(-1);
  }
  assert(mutexes.count(ml.addr));

  Mutex &mutex = mutexes[ml.addr];
  assert(0 <= mutex.last_access);
  assert(mutex.locked && "Unlocked mutex got unlocked again");

  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::M_UNLOCK;
  mayConflict(&ml);

  mutex.last_access = prefix_idx;
  mutex.locked = false;
}


void ZBuilderTSO::mutex_lock(const SymAddrSize &ml)
{
  //llvm::errs() << " M_LOCK_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << "_";
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)) {
    // Assume static initialization
    mutexes[ml.addr] = Mutex(-1);
  }
  assert(mutexes.count(ml.addr));
  Mutex &mutex = mutexes[ml.addr];

  assert(curnode().kind == ZEvent::Kind::DUMMY);
  // We make sure that every lock is an event with exactly
  // one instruction - only the lock instruction itself
  assert(prefix_idx == (int) (prefix.size() - 1));
  if (prefix[prefix_idx].size > 1) {
    // Some invisible instructions happened before the
    // lock instruction, we split this into two events
    unsigned p = curnode().iid.get_pid();
    --prefix[prefix_idx].size; // Subtract the lock instruction
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
    assert(prefix.back().traceID() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }

  assert(curnode().size == 1);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::M_LOCK;
  curnode().observed_trace_id = mutex.last_access; // initialized with -1

  mayConflict(&ml);
  if (!sch_replay)
    somethingToAnnotate.insert(curnode().iid.get_pid());
  endsWithLockFail.erase(curnode().iid.get_pid());

  assert(!mutex.locked);
  mutex.last_lock = mutex.last_access = prefix_idx;
  mutex.locked = true;
}


void ZBuilderTSO::mutex_lock_fail(const SymAddrSize &ml) {
  //llvm::errs() << " M_LOCKFAIL" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << "_";
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)){
    // Assume static initialization
    mutexes[ml.addr] = Mutex(-1);
  }
  assert(mutexes.count(ml.addr));
  #ifndef NDEBUG
  Mutex &mutex = mutexes[ml.addr];
  assert(0 <= mutex.last_lock);
  assert(mutex.locked);
  #endif

  if (!sch_replay)
    somethingToAnnotate.insert(curnode().iid.get_pid());
  endsWithLockFail.emplace(curnode().iid.get_pid(), ml);
}


std::pair<std::vector<ZEvent>, bool> ZBuilderTSO::extendGivenTrace() {
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


void ZBuilderTSO::add_failed_lock_attempts() {
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
    assert(prefix.back().traceID() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
    assert(curnode().size == 1);
    assert(curnode().kind == ZEvent::Kind::DUMMY);
    curnode().kind = ZEvent::Kind::M_LOCK;

    mayConflict(&(p_ml.second));
  }
}


Trace *ZBuilderTSO::get_trace() const
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
