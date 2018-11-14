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
  // Construct the explorer with:
  // the initial trace, this original TB, the star_root_index
  VCExplorer explorer = VCExplorer(std::move(prefix), *this, 1);
  assert(prefix.empty());

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

  AFTER_FAKING_LOCK:
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

    if (replay_trace[prefix_idx].kind == VCEvent::Kind::M_LOCKATTEMPT) {
      // Next follows a lock that is not annotated yet
      // We do not replay it, instead we push a lockattempt
      // event and stop this thread
      #ifndef NDEBUG
      // this is the last event of p
      for (unsigned i = prefix_idx+1; i < replay_trace.size(); ++i)
        assert(replay_trace[i].iid.get_pid() != (int) p);
      #endif
      prefix.emplace_back(IID<IPid>(IPid(p),threads[p].last_event_index()),
                          threads[p].cpid,
                          threads[p].executed_instructions + 1, // +1 so that first will be 1
                          threads[p].executed_events, // so that first will be 0
                          prefix.size(),
                          &replay_trace[prefix_idx]);
      prefix.back().setSize(replay_trace[prefix_idx].size);
      ++threads[p].executed_events;
      threads[p].executed_instructions += 1 + replay_trace[prefix_idx].size;
      assert(sch_replay);
      // stop this thread
      threads_with_unannotated_read.insert(curnode().iid.get_pid());
      assert(replay_trace[prefix_idx].size == prefix[prefix_idx].size);
      goto AFTER_FAKING_LOCK;
    }

    prefix.emplace_back(IID<IPid>(IPid(p),threads[p].last_event_index()),
                        threads[p].cpid,
                        threads[p].executed_instructions + 1, // +1 so that first will be 1
                        threads[p].executed_events, // so that first will be 0
                        prefix.size(),
                        nullptr);
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

  if (!in_critical_section.empty()) {
    assert(in_critical_section.size() == 1);
    unsigned cs_ipid = (unsigned) in_critical_section.begin()->first;
    if (!threads_with_unannotated_read.count(cs_ipid)) {
      bool ret = schedule_thread(proc, cs_ipid);
      assert(ret);
      return true;
    } else
      return false;
  }

  for(unsigned p = 0; p < threads.size(); p += 2) {
    if (!threads_with_unannotated_read.count(p)) {
      // We only schedule threads that have not yet seen a new visible read
      if (schedule_thread(proc, p))
        return true;
    }
  }

  // We did not schedule anything
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
                        prefix.size(),
                        nullptr);
    ++threads[p].executed_events;
    assert(prefix.back().id == prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
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
  if (ml)
    curn.ml = *ml;
  curn.instruction = current_inst;

  assert(!sch_replay ||
         replay_trace[prefix_idx].kind == prefix[prefix_idx].kind ||
         (replay_trace[prefix_idx].kind == VCEvent::Kind::M_LOCK &&
          prefix[prefix_idx].kind == VCEvent::Kind::M_LOCKATTEMPT));
}

void VCTraceBuilder::metadata(const llvm::MDNode *md)
{
  auto& cur = curnode();
  if (cur.md == nullptr)
    cur.md = md;
  last_md = md;
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
  // curnode().childs_cpid = CPS.dry_spawn(threads[parent_ipid].cpid);
  CPid child_cpid = CPS.spawn(threads[parent_ipid].cpid);
  curnode().childs_cpid = child_cpid;

  threads.emplace_back(child_cpid, prefix_idx); // second arg was threads[parent_ipid].event_indices
  threads.emplace_back(CPS.new_aux(child_cpid), prefix_idx); // second arg was threads[parent_ipid].event_indices
  threads.back().available = false; // Empty store buffer
}

void VCTraceBuilder::join(int tgt_proc)
{
  //llvm::errs() << " JOIN";
  assert(curnode().kind == VCEvent::Kind::DUMMY);
  curnode().kind = VCEvent::Kind::JOIN;
  mayConflict();
  curnode().childs_cpid = threads[2*tgt_proc].cpid;
}

void VCTraceBuilder::atomic_store(const SymData &sd, int val)
{
  // Stores to local memory on stack may not conflict
  assert(llvm::isa<llvm::StoreInst>(current_inst));
  //if (!llvm::isa<llvm::AllocaInst>(current_inst->getOperand(1)->stripInBoundsOffsets() )) {
  if (sd.get_ref().addr.block.is_global()) {
    assert(curnode().kind == VCEvent::Kind::DUMMY);
    curnode().kind = VCEvent::Kind::STORE;
    const SymAddrSize &ml = sd.get_ref();
    mayConflict(&ml);
    curnode().value = val;
    //llvm::errs() << " STORE_" << ml.to_string() << " ";
  }
}

void VCTraceBuilder::load(const SymAddrSize &ml, int val)
{
  // Loads from stack may not conflict
  assert(llvm::isa<llvm::LoadInst>(current_inst));
  //if (!llvm::isa<llvm::AllocaInst>(current_inst->getOperand(0)->stripInBoundsOffsets())) {
  if (ml.addr.block.is_global()) {
    assert(curnode().kind == VCEvent::Kind::DUMMY);
    curnode().kind = VCEvent::Kind::LOAD;
    mayConflict(&ml);
    curnode().value = val;
    //llvm::errs() << " LOAD_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << " ";
    if (sch_initial || sch_extend)
      if (threads.size() > 2)
        threads_with_unannotated_read.insert(curnode().iid.get_pid());
    // Otherwise the observation is determined and no need to annotate this node
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
  assert(mutexes.count(ml.addr) == 0);
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
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)){
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

  int p = curnode().iid.get_pid();
  assert(in_critical_section.count(p));
  if (in_critical_section[p] > 1)
      in_critical_section[p] =
        in_critical_section[p] - 1;
  else
    in_critical_section.erase(p);
}

void VCTraceBuilder::mutex_lock(const SymAddrSize &ml)
{
  //llvm::errs() << " M_LOCKATT_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << "_";
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)){
    // Assume static initialization
    mutexes[ml.addr] = Mutex();
  }
  assert(mutexes.count(ml.addr));
  Mutex &mutex = mutexes[ml.addr];

  assert(curnode().kind == VCEvent::Kind::DUMMY);
  curnode().kind = VCEvent::Kind::M_LOCKATTEMPT;

  mayConflict(&ml);
  curnode().value = mutex.value; // READ

  if (sch_initial || sch_extend) {
    // We see this for the first time,
    // note it down and stop here
    threads_with_unannotated_read.insert(curnode().iid.get_pid());
    doNotLock = true;
    //llvm::errs() << "NO ";
  }
  else {
    // This is a known and already annotated
    // lock, we proceed normally
    assert(!mutex.locked && mutex.value >= 0);
    mutex.last_lock = mutex.last_access = prefix_idx;
    mutex.locked = true;
    mutex.value = - 1 - (curnode().iid.get_pid() * 1000000);
                  // mutex locked by this process (<0)
    assert(sch_replay);
    assert(replay_trace[prefix_idx].kind == VCEvent::Kind::M_LOCK);
    curnode().kind = VCEvent::Kind::M_LOCK;
    int p = curnode().iid.get_pid();
    if (in_critical_section.count(p))
      in_critical_section[p] = in_critical_section[p] + 1;
    else
      in_critical_section.emplace(p, 1);
    doNotLock = false;
    //llvm::errs() << "YES ";
  }
}

void VCTraceBuilder::mutex_lock_fail(const SymAddrSize &ml){
  //llvm::errs() << " M_LOCKFAIL" << &ml;
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
}

std::pair<std::vector<VCEvent>,
          std::unordered_map<int, int>>
VCTraceBuilder::extendGivenTrace() {
  assert(sch_replay && !replay_trace.empty());

  std::unique_ptr<llvm::ExecutionEngine> EE(DPORDriver::create_execution_engine(M, *this, config));

  // Run main.
  EE->runFunctionAsMain(M->getFunction("main"), {"prog"}, 0);

  // Run static destructors.
  EE->runStaticConstructorsDestructors(true);

  return {prefix, in_critical_section};
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

  for(unsigned i = 0; i < prefix.size(); ++i){
    cmp.push_back(IID<CPid>(threads[prefix[i].iid.get_pid()].cpid,prefix[i].iid.get_index()));
    cmp_md.push_back(prefix[i].md);
  };
  for(unsigned i = 0; i < errors.size(); ++i){
    errs.push_back(errors[i]->clone());
  }
  Trace *t = new IIDSeqTrace(cmp,cmp_md,errs);
  t->set_blocked(false);
  return t;
}
