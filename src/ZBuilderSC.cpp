/* Copyright (C) 2016-2017 Marek Chalupa
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

#include "Debug.h"
#include "ZTrace.h"
#include "ZBuilderSC.h"
#include "ZExplorer.h"

static const bool DEBUG = false;
#include "ZDebug.h"


// Use at the very beginning to get an initial trace
ZBuilderSC::ZBuilderSC
(const Configuration &conf, llvm::Module *m, unsigned s_r_i)
  : TSOTraceBuilder(conf),
    sch_replay(false), sch_extend(true),
    star_root_index(s_r_i)
{
  config = &conf;
  M = m;
  prefix.reserve(64);
  ext_from_id = 0;
}


// Use when you want to do the following:
// (step1) replay the trace tr
// (step2) get a maximal extension
ZBuilderSC::ZBuilderSC
(const Configuration &conf, llvm::Module *m, std::vector<ZEvent>&& tr)
  : TSOTraceBuilder(conf),
    sch_replay(true), sch_extend(false),
    replay_trace(std::move(tr))
{
  config = &conf;
  M = m;
  prefix.reserve(replay_trace.size() + 64);
  assert(tr.empty() && !replay_trace.empty());
  ext_from_id = replay_trace.size();
}


bool ZBuilderSC::reset()
{
  #ifndef NDEBUG
    std::cout << "RUNNING DEBUG VERSION\n";
  #endif

  if (this->has_error()) {
    this->error_trace = this->get_trace();
    return true;
  }

  // Lock event for every thread ending with a failed mutex lock attempt
  add_failed_lock_attempts();

  // Construct the explorer with the original ZBuilder pointer
  ZExplorer explorer(*this);

  // Construct initial trace
  ZTrace initial_trace(
    nullptr, std::vector<ZEvent>(), std::vector<std::unique_ptr<ZEvent>>(),
    ZAnnotation(), std::set<ZEventID>(), ZEventID(true));

  // Construct initial extension
  ZTraceExtension initial_extension(
    std::move(prefix), ext_from_id,
    someThreadAssumeBlocked, !endsWithLockFail.empty());
  assert(prefix.empty());

  // Call the main method
  bool error = explorer.extend_and_explore(
    initial_trace, std::move(initial_extension));

  // Print the result statistics
  explorer.print_stats();

  return error;
}


/* *************************** */
/* SCHEDULING                  */
/* *************************** */


bool ZBuilderSC::schedule(int *proc, int *aux, int *alt, bool *dryrun)
{
  start_err("schedule");
  // Not using these arguments
  *dryrun = false;
  *alt = 0;
  this->dryrun = false;

  if(sch_replay) {
    assert(!sch_extend);
    auto res = schedule_replay_trace(proc, aux);
    end_err("?a");
    return res;
  }

  assert(!sch_replay && sch_extend);
  auto res = schedule_arbitrarily(proc, aux);
  end_err("?b");
  return res;
}


bool ZBuilderSC::schedule_replay_trace(int *proc, int *aux)
{
  start_err("schedule_replay_trace...");
  assert(!replay_trace.empty());
  assert(prefix_idx < (int) replay_trace.size());

  // prefix_idx is the index into prefix.
  // replay_trace_idx is the index into replay_trace
  // pointing to the event we are currently replaying.
  // They may not be equal (replay_trace_idx might be higher)
  // since buffer-writes in replay_trace are skipped and
  // replayed later in order to directly precede memory-writes

  err_msg("prefix = " + std::to_string(prefix_idx) + " vs. replay_trace_size = " + std::to_string(replay_trace.size()));
  if (prefix_idx == -1 || replay_trace[replay_trace_idx].size == prefix[prefix_idx].size) {
    ++prefix_idx;
    ++replay_trace_idx;
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
      auto res = schedule_arbitrarily(proc, aux);
      end_err("?a");
      return res;
    }
    // Next instruction is the beginning of a new event
    assert(replay_trace.size() > prefix.size());
    assert(prefix_idx == (int) prefix.size());
    // Skip buffer-write events of replay_trace; each will be replayed
    // later in order to directly precede its memory-write counterpart
    unsigned p = replay_trace[replay_trace_idx].iid.get_pid();
    while (isWriteB(replay_trace[replay_trace_idx])) {
      assert(replay_trace[replay_trace_idx].instruction_id ==
             threads[p].executed_instructions + 1 && "Inconsistent scheduling");
      assert(replay_trace[replay_trace_idx].event_id() ==
             threads[p].executed_events && "Inconsistent scheduling");
      assert(!delayed_bwrite.count(p));
      delayed_bwrite.emplace(p);
      ++replay_trace_idx;
      assert(replay_trace_idx < replay_trace.size());
      p = replay_trace[replay_trace_idx].iid.get_pid();
    }
    // If below assertion fails, search for p' tied to r_t[p_i].cpid,
    // scan r_t down and swap p with p' in all found events
    assert(replay_trace[replay_trace_idx].cpid() ==
           threads[p].cpid && "IPID<->CPID correspondence has changed");
    if (isWriteM(replay_trace[replay_trace_idx])) {
      // We reached a memory-write in replay_trace, so now we will
      // replay the main thread to perform the buffer-write part
      // (delayed until now). Set ipid to the main-tread one; the size
      // in the memory-write event matches the desired size of buffer-write,
      // so replay that many steps; in the last replay step of buffer-write,
      // right after it the memory-write event gets created too.
      assert(p % 2 == 1);
      --p;
      assert(delayed_bwrite.count(p));
      delayed_bwrite.erase(p);
      // Mark that thread p owns the new event
      threads[p].event_indices.push_back(prefix_idx);
    } else {
      // Mark that thread p owns the new event
      threads[p].event_indices.push_back(prefix_idx);
      assert(replay_trace[replay_trace_idx].instruction_id ==
             threads[p].executed_instructions + 1 && "Inconsistent scheduling");
      assert(replay_trace[replay_trace_idx].event_id() ==
             threads[p].executed_events && "Inconsistent scheduling");
    }
    assert(!delayed_bwrite.count(p));
    // Create the new event
    prefix.emplace_back(IID<IPid>(IPid(p),threads[p].last_event_index()),
                        threads[p].cpid,
                        threads[p].executed_instructions + 1, // +1 so that first will be 1
                        threads[p].executed_events, // so that first will be 0
                        prefix.size());
    curnode_set_thread_id();
    fence(); // Each event in SC has fence
    // Mark that thread executes a new event
    ++threads[p].executed_events;
    assert(prefix.back().trace_id() == (int) prefix.size() - 1);
  } else {
    // Next instruction is a continuation of the current event
    assert(replay_trace[replay_trace_idx].size > prefix[prefix_idx].size);
    // Increase the size of the current event, it will be scheduled
    ++prefix[prefix_idx].size;
    assert(replay_trace[replay_trace_idx].size >=
           prefix[prefix_idx].size && "Inconsistent scheduling");
  }

  assert((unsigned) prefix_idx < replay_trace.size());
  unsigned p = replay_trace[replay_trace_idx].iid.get_pid();
  if (isWriteM(replay_trace[replay_trace_idx])) {
    // Set ipid to the main-tread one
    assert(p % 2 == 1);
    --p;
  }
  err_msg("Thread " + std::to_string(p) + ", instruction " + std::to_string(threads[p].executed_instructions));
  // Mark that thread p executes a new instruction
  ++threads[p].executed_instructions;
  assert(threads[p].available);

  // Here used to be:
  // ++threads[p].clock[p];
  // After refactoring, it should be:
  // threads[p].event_indices.push_back(prefix_idx);
  // But I don't know why here, I would suspect to
  // run this only when we have a new event

  bool ret = schedule_thread(proc, aux, p);
  assert(ret && "Bug in scheduling: could not reproduce a given replay trace");

  end_err("?b");
  return ret;
}


bool ZBuilderSC::schedule_arbitrarily(int *proc, int *aux)
{
  start_err("schedule_arbitrarily...");
  assert(!sch_replay && sch_extend);

  assert(threads.size() % 2 == 0);
  // We prefer auxiliary before real threads

  const unsigned sz = threads.size();
  unsigned p;

  #ifndef NDEBUG
  for (p = 1; p < sz; p += 2) { // Loop through auxiliary threads
    assert(!threads[p].available);
  }
  #endif

  for (p = 0; p < sz; p += 2) { // Loop through real threads
    if (schedule_thread(proc, aux, p)) {
      end_err("1d");
      return true;
    }
  }

  // We did not schedule anything
  prefix.shrink_to_fit();
  end_err("0");
  return false;
}


// Schedule the next instruction to be the next instruction from thread p
bool ZBuilderSC::schedule_thread(int *proc, int *aux, unsigned p)
{
  start_err("schedule_thread...");
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
    assert(*aux == -1 && "Only main threads are played in SC");

    end_err("1");
    return true;
  }

  end_err("0");
  return false;
}


// Used only when we are not replaying replay_trace
void ZBuilderSC::update_prefix(unsigned p)
{
  // Here used to be:
  // ++threads[p].clock[p];
  // After refactoring, it should be:
  // threads[p].event_indices.push_back(prefix_idx);
  // But I don't know why here, I would suspect to
  // run this only when we have a new event
  assert(sch_extend && !sch_replay);
  assert(p % 2 == 0);
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
    curnode_set_thread_id();
    ++threads[p].executed_events;
    fence(); // Each event in SC has fence
    assert(prefix.back().trace_id() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }
}


void ZBuilderSC::curnode_set_thread_id()
{
  unsigned thr_idx = INT_MAX;

  // Check if this process is already known
  auto ipidit = ipid_to_thraux.find(curnode().iid.get_pid());
  if (ipidit != ipid_to_thraux.end()) {
    thr_idx = ipidit->second.first;
    assert(curnode().aux_id() == ipidit->second.second);
  } else {
    auto it = proc_seq_to_thread_id.find(curnode().cpid().get_proc_seq());
    if (it == proc_seq_to_thread_id.end()) {
      thr_idx = proc_seq_to_thread_id.size();
      proc_seq_to_thread_id.emplace_hint(
        it, curnode().cpid().get_proc_seq(), thr_idx);
    } else
      thr_idx = it->second;

    // add to ipid cache for faster lookup next time
    ipid_to_thraux.emplace(
      curnode().iid.get_pid(),
      std::pair<unsigned, int>(thr_idx, curnode().aux_id()));
  }

  assert(thr_idx < 10000);
  curnode()._thread_id = (int) thr_idx;
}


/* ******************************************************* */
/* ADDRESSING FEEDBACK FROM INTERPRETER AFTER SCHEDULING   */
/* ******************************************************* */


void ZBuilderSC::refuse_schedule()
{
  start_err("refuse_schedule...");
  //llvm::errs() << " REFUSESCH";
  assert(prefix_idx == int(prefix.size())-1);
  assert(!prefix.back().may_conflict);
  assert(sch_extend && !sch_replay &&
         "Refusing during replay means inconsistent replay trace");

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
  end_err();
}


void ZBuilderSC::mayConflict(const SymAddrSize *ml)
{
  auto& curn = curnode();
  curn.may_conflict = true;
  if (ml) curn.ml = *ml;

  #ifndef NDEBUG
  bool consistent = (!sch_replay ||
                     replay_trace[replay_trace_idx].kind == prefix[prefix_idx].kind ||
                     (isWriteM(replay_trace[replay_trace_idx]) &&
                      isWriteB(prefix[prefix_idx])));
  if (!consistent) {
    llvm::errs() << "TRACE_TO_REPLAY\n";
    dump_trace(replay_trace);
    llvm::errs() << "\nACTUALLY_REPLAYED\n";
    dump_trace(prefix);
    llvm::errs() << "Problematic event (didn't go any further):\n";
    llvm::errs() << "TRACE_TO_REPLAY[" << prefix_idx << "]   ::: ";
    replay_trace[replay_trace_idx].dump();
    llvm::errs() << "ACTUALLY_REPLAYED[" << prefix_idx << "] ::: ";
    prefix[prefix_idx].dump();
  }
  assert(consistent);
  #endif
}


void ZBuilderSC::metadata(const llvm::MDNode *md)
{
  if (md) curnode().md = md;
}


IID<CPid> ZBuilderSC::get_iid() const
{
  IPid pid = curnode().iid.get_pid();
  int idx = curnode().iid.get_index();
  return IID<CPid>(threads[pid].cpid,idx);
}


void ZBuilderSC::spawn()
{
  start_err("spawn");
  //llvm::errs() << " SPAWN";
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::SPAWN;
  // store the CPid of the new thread
  IPid parent_ipid = curnode().iid.get_pid();
  CPid child_cpid = CPS.spawn(threads[parent_ipid].cpid);
  curnode().childs_cpid = child_cpid;
  mayConflict();

  threads.emplace_back(child_cpid, prefix_idx); // second arg was threads[parent_ipid].event_indices
  threads.emplace_back(CPS.new_aux(child_cpid), prefix_idx); // second arg was threads[parent_ipid].event_indices
  threads.back().available = false; // New thread starts with an empty store buffer
  end_err();
}


void ZBuilderSC::join(int tgt_proc)
{
  start_err("join...");
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
    assert(sch_extend && !sch_replay &&
           "Join with invisible instruction(s) should not be one event in replay trace");
    unsigned p = curnode().iid.get_pid();
    --prefix[prefix_idx].size; // Subtract the join instruction
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
    curnode_set_thread_id();
    fence(); // Each event in SC has fence
    ++threads[p].executed_events;
    assert(prefix.back().trace_id() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }

  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::JOIN;
  curnode().childs_cpid = threads[2*tgt_proc].cpid;
  mayConflict();
  end_err();
}


void ZBuilderSC::store(const SymData &sd, int val)
{
  start_err("store...");
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);

  unsigned p = curnode().iid.get_pid();
  assert(p % 2 == 0);

  // Do not put store to buffer in SC
  //threads[p].store_buffer.push_back(PendingStore(sd.get_ref(),prefix_idx,last_md));
  //threads[p+1].available = true;

  const SymAddrSize& ml = sd.get_ref();
  // Stores to local memory on stack may not conflict
  assert(ml.addr.block.is_global() || ml.addr.block.is_heap());

  // Global stores happening with just one thread existing may not conflict
  // but we have to record them as such anyway to keep track of them

  //llvm::errs() << " BUFFERWRITE_" << ml.to_string() << " ";
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::WRITEB;
  curnode().value = val;
  mayConflict(&ml);

  // Do not put store to visible store queue in SC
  // Create memory-write counterpart
  atomic_store(sd);
  end_err();
}


void ZBuilderSC::atomic_store(const SymData &sd)
{
  start_err("store...");
  assert(!dryrun);

  assert(isWriteB(curnode()));
  unsigned realp = curnode().iid.get_pid();
  assert(realp % 2 == 0 && "This function got called from buffer-write counterpart");
  unsigned auxp = realp + 1;
  int buffer_size = curnode().size; // We copy the size of buffer-write counterpart

  const SymAddrSize& ml = sd.get_ref();
  // We consider all Stack writes as invisible and atomic,
  // hence we do not put them into the store buffer
  assert((ml.addr.block.is_global() || ml.addr.block.is_heap()) &&
         "No Stack writes are allowed in the store buffer");

  // Store buffer is always empty in SC
  assert(threads[realp].store_buffer.empty() && !threads[auxp].available);

  // Create a new event for the memory-write counterpart
  ++prefix_idx;
  // Mark that auxiliary thread auxp owns the new event
  threads[auxp].event_indices.push_back(prefix_idx);
  // Create the new event
  prefix.emplace_back(IID<IPid>(IPid(auxp),threads[auxp].last_event_index()),
                      threads[auxp].cpid,
                      threads[auxp].executed_instructions, // first will be 1
                      threads[auxp].executed_events, // first will be 0
                      prefix.size());
  curnode_set_thread_id();
  // Auxiliary events do not have fence themselves -- fence();
  ++threads[auxp].executed_events;
  curnode().size = buffer_size; // We copy the size of buffer-write counterpart
  assert(prefix.back().trace_id() == (int) prefix.size() - 1);
  assert((unsigned) prefix_idx == prefix.size() - 1);

  //llvm::errs() << " MEMORYWRITE_" << ml.to_string() << " ";
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::WRITEM;
  mayConflict(&ml);

  // Corresponding buffer-write is the preceding event in SC
  curnode().write_other_trace_id = prefix_idx - 1;
  assert(curnode().write_other_trace_id >= 0 &&
         curnode().write_other_trace_id == (int) prefix.size() - 2);
  assert(prefix[curnode().write_other_trace_id].kind == ZEvent::Kind::WRITEB);
  assert(prefix[curnode().write_other_trace_id].ml == ml &&
         "Inconsistent memory locations");
  curnode().value = prefix[curnode().write_other_trace_id].value;
  assert(prefix[curnode().write_other_trace_id].write_other_trace_id == -1);
  prefix[curnode().write_other_trace_id].write_other_trace_id = prefix_idx;

  // Visible store queue is always empty in SC
  assert(visibleStoreQueue.empty());

  lastWrite[ml] = prefix_idx;
  end_err();
}


void ZBuilderSC::load(const SymAddrSize &ml, int val)
{
  start_err("load");
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
    // no buffer-write queue in SC
    assert(visibleStoreQueue.empty());
    curnode().observed_trace_id = obs_idx;
  }
  end_err();
}


void ZBuilderSC::fence()
{
  start_err("fence");
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  assert(!visibleStoreQueue.count(curnode().iid.get_pid()) ||
         visibleStoreQueue[curnode().iid.get_pid()].empty());
  assert(threads[curnode().iid.get_pid()].store_buffer.empty());
  curnode().fence = true;
  end_err();
}


void ZBuilderSC::mutex_init(const SymAddrSize &ml)
{
  start_err("mutex_init");
  //llvm::errs() << " M_INIT_" << ml.to_string() << " ";
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  fence();
  assert(!mutexes.count(ml.addr));
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::M_INIT;
  mayConflict(&ml);
  mutexes[ml.addr] = Mutex(-1); // prefix_idx
  end_err();
}


void ZBuilderSC::mutex_destroy(const SymAddrSize &ml)
{
  start_err("mutex_destroy");
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
  end_err();
}


void ZBuilderSC::mutex_unlock(const SymAddrSize &ml)
{
  start_err("mutex_unlock");
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
  end_err();
}


void ZBuilderSC::mutex_lock(const SymAddrSize &ml)
{
  start_err("mutex_lock");
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
    assert(sch_extend && !sch_replay &&
           "Lock with invisible instruction(s) should not be one event in replay trace");
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
    curnode_set_thread_id();
    fence(); // Each event in SC has fence
    ++threads[p].executed_events;
    assert(prefix.back().trace_id() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }

  assert(curnode().size == 1);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::M_LOCK;
  curnode().observed_trace_id = mutex.last_access; // initialized with -1

  mayConflict(&ml);
  endsWithLockFail.erase(curnode().iid.get_pid());

  assert(!mutex.locked);
  mutex.last_lock = mutex.last_access = prefix_idx;
  mutex.locked = true;
  end_err();
}


void ZBuilderSC::mutex_lock_fail(const SymAddrSize &ml) {
  start_err("mutex_lock_fail");
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

  endsWithLockFail.emplace(curnode().iid.get_pid(), ml);
  end_err();
}


ZTraceExtension ZBuilderSC::extendGivenTrace() {
  start_err("extendGivenTrace...");
  assert(sch_replay && !replay_trace.empty());

  std::unique_ptr<llvm::ExecutionEngine> EE(DPORDriver::create_execution_engine(M, *this, *config));

  // Run main.
  EE->runFunctionAsMain(M->getFunction("main"), {"prog"}, 0);

  // Run static destructors.
  EE->runStaticConstructorsDestructors(true);

  // Lock event for every thread ending with a failed mutex lock attempt
  add_failed_lock_attempts();

  if (has_error()) {
    end_err("found-error");
    return ZTraceExtension();
  }

  ZTraceExtension res(
    std::move(prefix), ext_from_id,
    someThreadAssumeBlocked, !endsWithLockFail.empty());
  assert(prefix.empty());
  end_err();
  return res;
}


void ZBuilderSC::add_failed_lock_attempts() {
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
    curnode_set_thread_id();
    ++threads[p].executed_events;
    assert(prefix.back().trace_id() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
    assert(curnode().size == 1);
    assert(curnode().kind == ZEvent::Kind::DUMMY);
    curnode().kind = ZEvent::Kind::M_LOCK;
    curnode().failed_lock = true;

    fence(); // Each event in SC has fence
    mayConflict(&(p_ml.second));
  }
}


Trace *ZBuilderSC::get_trace() const
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
