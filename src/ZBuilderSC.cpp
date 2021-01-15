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
(const Configuration &conf, llvm::Module *m, bool init_only)
  : TSOTraceBuilder(conf),
    sch_replay(false), sch_extend(true),
    prefix(new std::vector<std::unique_ptr<ZEvent>>()),
    initial_trace_only(init_only),
    graph(new ZGraph()),
    po_part(new ZPartialOrder(*graph)),
    po_full(new ZPartialOrder(*graph))
{
  config = &conf;
  M = m;
  prefix->reserve(64);
}


// Use when you want to do the following:
// (step1) replay the trace tr
// (step2) get a maximal extension
ZBuilderSC::ZBuilderSC
(const Configuration &conf, llvm::Module *m, std::vector<ZEvent>&& tr,
 ZPartialOrder&& mutated_po)
  : TSOTraceBuilder(conf),
    sch_replay(true), sch_extend(false),
    prefix(new std::vector<std::unique_ptr<ZEvent>>()),
    replay_trace(std::move(tr)),
    graph(new ZGraph()),
    po_part(new ZPartialOrder(std::move(mutated_po), *graph)),
    po_full(new ZPartialOrder(*graph))
{
  config = &conf;
  M = m;
  prefix->reserve(replay_trace.size() + 64);
}


bool ZBuilderSC::reset()
{
  if (this->has_error()) {
    this->error_trace = this->get_trace();
    if (initial_trace_only) dump_trace(*prefix);
    return true;
  }

  if (initial_trace_only) {
    dump_trace(*prefix);
    return false;
  }

  // Add lock event for every thread ending with a failed mutex lock attempt
  add_failed_lock_attempts();

  // Graph and PO building - process the very last event
  graph_po_process_event(true);
  assert(graph_po_constructed);

  // Construct the explorer with the original TB, holding the initial trace
  ZExplorer explorer = ZExplorer(*this);

  // Call the main method
  bool error = explorer.explore();

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
  *aux = -1;
  this->dryrun = false;

  if(sch_replay) {
    assert(!sch_extend);
    auto res = schedule_replay_trace(proc);
    end_err("?a");
    return res;
  }

  assert(!sch_replay && sch_extend);
  auto res = schedule_arbitrarily(proc);
  end_err("?b");
  return res;
}


bool ZBuilderSC::schedule_replay_trace(int *proc)
{
  start_err("schedule_replay_trace...");
  assert(!replay_trace.empty());
  assert(prefix_idx < (int) replay_trace.size());

  // prefix_idx is the index into prefix. Since we are up until
  // now only replaying the event sequence of replay_trace,
  // prefix_idx is also the index into replay_trace pointing
  // to the event we are currently replaying

  err_msg("prefix = " + std::to_string(prefix_idx) + " vs. replay_trace_size = " + std::to_string(replay_trace.size()));
  if (prefix_idx == -1 || replay_trace[prefix_idx].size == (*prefix)[prefix_idx]->size) {
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
      auto res = schedule_arbitrarily(proc);
      end_err("?a");
      return res;
    }
    // Next instruction is the beginning of a new event
    assert(replay_trace.size() > prefix->size());
    assert(prefix_idx == (int) prefix->size());
    unsigned p = replay_trace[prefix_idx].iid.get_pid();
    // If below assertion fails, search for p' tied to r_t[p_i].cpid,
    // scan r_t down and swap p with p' in all found events
    assert(replay_trace[prefix_idx].cpid() ==
           threads[p].cpid && "IPID<->CPID correspondence has changed");
    // Mark that thread p owns the new event
    threads[p].event_indices.push_back(prefix_idx);
    assert(replay_trace[prefix_idx].instruction_id ==
           threads[p].executed_instructions + 1 && "Inconsistent scheduling");
    assert(replay_trace[prefix_idx].event_id() ==
           threads[p].executed_events && "Inconsistent scheduling");
    // Create the new event
    prefix->emplace_back(new ZEvent(
      IID<IPid>(IPid(p),threads[p].last_event_index()),
      threads[p].cpid,
      threads[p].executed_instructions + 1, // +1 so that first will be 1
      threads[p].executed_events, // so that first will be 0
      prefix->size()));
    // Mark that thread executes a new event
    ++threads[p].executed_events;
    assert(prefix->back()->trace_id() == (int) prefix->size() - 1);
    // Graph and PO building
    if (prefix->size() != 1) graph_po_process_event(false);
  } else {
    // Next instruction is a continuation of the current event
    assert(replay_trace[prefix_idx].size > (*prefix)[prefix_idx]->size);
    // Increase the size of the current event, it will be scheduled
    ++((*prefix)[prefix_idx]->size);
    assert(replay_trace[prefix_idx].size >=
           (*prefix)[prefix_idx]->size && "Inconsistent scheduling");
  }

  assert((unsigned) prefix_idx < replay_trace.size());
  unsigned p = replay_trace[prefix_idx].iid.get_pid();
  err_msg("Thread " + std::to_string(p) + ", instruction " + std::to_string(threads[p].executed_instructions));
  // Mark that thread p executes a new instruction
  ++threads[p].executed_instructions;
  assert(threads[p].available);

  bool ret = schedule_thread(proc, p);
  assert(ret && "Bug in scheduling: could not reproduce a given replay trace");

  end_err("?b");
  return ret;
}


bool ZBuilderSC::schedule_arbitrarily(int *proc)
{
  start_err("schedule_arbitrarily...");
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
    #ifndef NDEBUG
    for (p = 1; p < sz; p += 2) { // Loop through auxiliary threads
      assert(!threads[p].available);
    }
    #endif
    for (p = 0; p < sz; p += 2) { // Loop through real threads
      if (!somethingToAnnotate.count(p))
        if (schedule_thread(proc, p)) {
          end_err("1b");
          return true;
        }
    }
  }

  #ifndef NDEBUG
  for (p = 1; p < sz; p += 2) { // Loop through auxiliary threads
    assert(!threads[p].available);
  }
  #endif

  for (p = 0; p < sz; p += 2) { // Loop through real threads
    if (schedule_thread(proc, p)) {
      end_err("1d");
      return true;
    }
  }

  // We did not schedule anything
  prefix->shrink_to_fit();
  end_err("0");
  return false;
}


// Schedule the next instruction to be the next instruction from thread p
bool ZBuilderSC::schedule_thread(int *proc, unsigned p)
{
  start_err("schedule_thread...");
  if (threads[p].available && !threads[p].sleeping &&
      (conf.max_search_depth < 0 || threads[p].last_event_index() < conf.max_search_depth)) {

    if (sch_extend) {
      update_prefix(p);
    } else {
      assert(sch_replay);
    }

    // set the scheduled thread
    *proc = p/2;

    end_err("1");
    return true;
  }

  end_err("0");
  return false;
}


// Used only when we are not replaying replay_trace
void ZBuilderSC::update_prefix(unsigned p)
{
  start_err("update_prefix...");
  assert(sch_extend && !sch_replay);
  assert(p % 2 == 0);
  assert(p < threads.size());
  ++threads[p].executed_instructions;

  if (prefix_idx != -1 && (int) p == curnode().iid.get_pid() &&
      !curnode().may_conflict) {
    // 1) This will not be the very first instruction
    // 2) Context switch is not going to happen now
    // 3) In this event we have only invisible instructions so far
    // Because of the above, we extend the current event
    assert(prefix_idx == (int) (prefix->size() - 1));
    ++((*prefix)[prefix_idx]->size);
    end_err("event_extend");
  } else {
    // Because one of 1)2)3) doesn't hold, we create a new event
    ++prefix_idx;
    // Mark that thread p owns the new event
    threads[p].event_indices.push_back(prefix_idx);
    // Create the new event
    prefix->emplace_back(new ZEvent(
      IID<IPid>(IPid(p),threads[p].last_event_index()),
      threads[p].cpid,
      threads[p].executed_instructions, // first will be 1
      threads[p].executed_events, // first will be 0
      prefix->size()));
    ++threads[p].executed_events;
    assert(prefix->back()->trace_id() == (int) prefix->size() - 1);
    assert((unsigned) prefix_idx == prefix->size() - 1);
    // Graph and PO building
    if (prefix->size() != 1) graph_po_process_event(false);
    end_err("event_new");
  }
}


/* ******************************************************* */
/* ADDRESSING FEEDBACK FROM INTERPRETER AFTER SCHEDULING   */
/* ******************************************************* */


void ZBuilderSC::refuse_schedule()
{
  start_err("refuse_schedule...");
  assert(prefix_idx == int(prefix->size())-1);
  assert(!prefix->back()->may_conflict);
  assert(sch_extend && !sch_replay &&
         "Refusing during replay means inconsistent replay trace");

  unsigned p = curnode().iid.get_pid();
  assert(p % 2 == 0 && "Only real threads can be refused");

  --threads[p].executed_instructions; //

  if (curnode().size == 1) {
    // Refused instruction wanted to create a new event,
    // therefore this entire event needs to be deleted
    // Unmark ownership of this new event
    assert((int) threads[p].event_indices.back() == prefix_idx);
    threads[p].event_indices.pop_back();
    --threads[p].executed_events; //
    // Delete this new event
    prefix->pop_back();
    --prefix_idx;
  } else {
    // Refused instruction wanted to extend the current event
    --curnode().size;
  }

  assert(prefix_idx == int(prefix->size())-1);

  // Mark this thread as unavailable since its next instruction
  // is the one we tried to schedule now, and it got refused
  mark_unavailable(p/2, -1);
  end_err();
}


void ZBuilderSC::mayConflict(const SymAddrSize *ml)
{
  auto& curn = curnode();
  curn.may_conflict = true;
  if (ml) curn._ml = *ml;

  #ifndef NDEBUG
  bool consistent = (!sch_replay ||
                     replay_trace[prefix_idx].kind == (*prefix)[prefix_idx]->kind);
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
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::SPAWN;
  // store the CPid of the new thread
  IPid parent_ipid = curnode().iid.get_pid();
  CPid child_cpid = CPS.spawn(threads[parent_ipid].cpid);
  curnode()._childs_cpid = child_cpid;
  mayConflict();

  threads.emplace_back(child_cpid, prefix_idx); // second arg was threads[parent_ipid].event_indices
  threads.emplace_back(CPS.new_aux(child_cpid), prefix_idx); // second arg was threads[parent_ipid].event_indices
  threads.back().available = false; // New thread starts with an empty store buffer
  end_err();
}


void ZBuilderSC::join(int tgt_proc)
{
  start_err("join...");
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  // We make sure that every join is an event with exactly
  // one instruction - only the join instruction itself
  assert(prefix_idx == (int) (prefix->size() - 1));
  if ((*prefix)[prefix_idx]->size > 1) {
    // Some invisible instructions happened before the
    // join instruction, we split this into two events
    assert(sch_extend && !sch_replay &&
           "Join with invisible instruction(s) should not be one event in replay trace");
    unsigned p = curnode().iid.get_pid();
    --((*prefix)[prefix_idx]->size); // Subtract the join instruction
    // Create a new event with only the join instruction
    ++prefix_idx;
    // Mark that thread p owns the new event
    threads[p].event_indices.push_back(prefix_idx);
    // Create the new event
    prefix->emplace_back(new ZEvent(
      IID<IPid>(IPid(p),threads[p].last_event_index()),
      threads[p].cpid,
      threads[p].executed_instructions, // first will be 1
      threads[p].executed_events, // first will be 0
      prefix->size()));
    ++threads[p].executed_events;
    assert(prefix->back()->trace_id() == (int) prefix->size() - 1);
    assert((unsigned) prefix_idx == prefix->size() - 1);
    // Graph and PO building
    graph_po_process_event(false);
  }

  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::JOIN;
  curnode()._childs_cpid = threads[2*tgt_proc].cpid;
  mayConflict();
  end_err();
}


void ZBuilderSC::atomic_store(const SymAddrSize &ml, int val)
{
  start_err("store...");
  assert(!dryrun);

  unsigned realp = curnode().iid.get_pid();
  assert(realp % 2 == 0);
  assert(threads[realp].store_buffer.empty() && !threads[realp + 1].available);

  // We consider all Stack writes as invisible
  assert((ml.addr.block.is_global() || ml.addr.block.is_heap()) &&
         "Stack writes are treated as invisible");

  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::WRITE;
  curnode()._value = val;
  mayConflict(&ml);

  assert(prefix->back()->trace_id() == (int) prefix->size() - 1);
  assert(prefix_idx == (int) prefix->size() - 1);

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
      threads.size() > 2) { // 2 because there are two entries for each thread
    assert(curnode().kind == ZEvent::Kind::DUMMY);
    curnode().kind = ZEvent::Kind::READ;
    curnode()._value = val;
    mayConflict(&ml);

    // Set ID of observed event
    int obs_idx = (lastWrite.count(ml)) ? lastWrite[ml] : -1;
    curnode()._observed_trace_id = obs_idx;

    if (!sch_replay) {
      somethingToAnnotate.insert(curnode().iid.get_pid());
      threads_with_unannotated_readlock.insert(curnode().cpid());
    }
  }
  end_err();
}


void ZBuilderSC::compare_exchange
(const SymAddrSize &ml, int old_val, int compare_val, int exchange_val, bool success)
{
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  assert(success || old_val != compare_val);
  assert(!success || old_val == compare_val);
  if (!ml.addr.block.is_global() && !ml.addr.block.is_heap())
    return;
  if (threads.size() <= 2) { // 2 because there are two entries for each thread
    // Only record the write (in case it happened)
    if (success) {
      assert(curnode().kind == ZEvent::Kind::DUMMY);
      atomic_store(ml, exchange_val);
      assert(is_write(curnode()));
    }
    return;
  }
  start_err("compare_exchange");

  if (sch_replay) {
    assert(prefix_idx < (int) replay_trace.size());
    assert(replay_trace.at(prefix_idx).is_read_of_cas());
    assert(replay_trace.at(prefix_idx).size > 0);
    if (success) {
      assert(prefix_idx + 1 < (int) replay_trace.size());
      assert(is_write(replay_trace.at(prefix_idx + 1)));
      assert(same_ml(replay_trace.at(prefix_idx), replay_trace.at(prefix_idx + 1)));
      assert(replay_trace.at(prefix_idx).cpid() == replay_trace.at(prefix_idx + 1).cpid());
      assert(replay_trace.at(prefix_idx + 1).is_write_of_cas());
      assert(replay_trace.at(prefix_idx + 1).size == 0);
    } else {
      assert(prefix_idx + 1 >= (int) replay_trace.size() ||
             !replay_trace.at(prefix_idx + 1).is_write_of_cas());
      assert(prefix_idx + 1 >= (int) replay_trace.size() ||
             replay_trace.at(prefix_idx + 1).size > 0);
    }
  }

  load(ml, old_val);
  curnode().set_read_of_cas(compare_val, exchange_val);

  if (!success) {
    end_err();
    return;
  }

  // Create a write event of instruction size 0,
  // it is the continuation of the atomic CAS
  assert(prefix_idx == (int) (prefix->size() - 1));
  unsigned p = curnode().iid.get_pid();
  ++prefix_idx;
  // Mark that thread p owns the new event
  threads[p].event_indices.push_back(prefix_idx);
  // Create the new event
  prefix->emplace_back(new ZEvent(
    IID<IPid>(IPid(p),threads[p].last_event_index()),
    threads[p].cpid,
    threads[p].executed_instructions, // first will be 1
    threads[p].executed_events, // first will be 0
    prefix->size()));
  ++threads[p].executed_events;
  assert(prefix->back()->trace_id() == (int) prefix->size() - 1);
  assert((unsigned) prefix_idx == prefix->size() - 1);
  assert(curnode().kind == ZEvent::Kind::DUMMY);

  atomic_store(ml, exchange_val);
  assert(is_write(curnode()));
  curnode().size = 0;
  curnode().set_write_of_cas();
  // Graph and PO building
  graph_po_process_event(false);

  end_err();
}


void ZBuilderSC::read_modify_write
(const SymAddrSize &ml, int old_val, int new_val)
{
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  if (!ml.addr.block.is_global() && !ml.addr.block.is_heap())
    return;
  if (threads.size() <= 2) { // 2 because there are two entries for each thread
    // Only record the write
    assert(curnode().kind == ZEvent::Kind::DUMMY);
    atomic_store(ml, new_val);
    assert(is_write(curnode()));
    return;
  }
  if (threads.size() <= 2)
    return; // 2 because there are two entries for each thread
  start_err("read_modify_write");

  if (sch_replay) {
    assert(prefix_idx + 1 < (int) replay_trace.size());
    assert(replay_trace.at(prefix_idx).is_read_of_rmw());
    assert(replay_trace.at(prefix_idx).size > 0);
    assert(is_write(replay_trace.at(prefix_idx + 1)));
    assert(same_ml(replay_trace.at(prefix_idx), replay_trace.at(prefix_idx + 1)));
    assert(replay_trace.at(prefix_idx).cpid() == replay_trace.at(prefix_idx + 1).cpid());
    assert(replay_trace.at(prefix_idx + 1).is_write_of_rmw());
    assert(replay_trace.at(prefix_idx + 1).size == 0);
  }

  load(ml, old_val);
  curnode().set_read_of_rmw();

  // Create a write event of instruction size 0,
  // it is the continuation of the atomic RMW
  assert(prefix_idx == (int) (prefix->size() - 1));
  unsigned p = curnode().iid.get_pid();
  ++prefix_idx;
  // Mark that thread p owns the new event
  threads[p].event_indices.push_back(prefix_idx);
  // Create the new event
  prefix->emplace_back(new ZEvent(
    IID<IPid>(IPid(p),threads[p].last_event_index()),
    threads[p].cpid,
    threads[p].executed_instructions, // first will be 1
    threads[p].executed_events, // first will be 0
    prefix->size()));
  ++threads[p].executed_events;
  assert(prefix->back()->trace_id() == (int) prefix->size() - 1);
  assert((unsigned) prefix_idx == prefix->size() - 1);
  assert(curnode().kind == ZEvent::Kind::DUMMY);

  atomic_store(ml, new_val);
  assert(is_write(curnode()));
  curnode().size = 0;
  curnode().set_write_of_rmw();
  // Graph and PO building
  graph_po_process_event(false);

  end_err();
}


void ZBuilderSC::fence()
{
  //start_err("fence");
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  assert(threads[curnode().iid.get_pid()].store_buffer.empty());
  //end_err();
}


void ZBuilderSC::mutex_init(const SymAddrSize &ml)
{
  start_err("mutex_init");
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  fence();
  assert(!mutexes.count(ml.addr));
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  // curnode().kind = ZEvent::Kind::M_INIT;
  // mayConflict(&ml);
  mutexes[ml.addr] = Mutex(-1); // prefix_idx
  end_err();
}


void ZBuilderSC::mutex_destroy(const SymAddrSize &ml)
{
  start_err("mutex_destroy");
  assert(!dryrun);
  assert(curnode().iid.get_pid() % 2 == 0);
  fence();
  if(!conf.mutex_require_init && !mutexes.count(ml.addr)) {
    // Assume static initialization
    mutexes[ml.addr] = Mutex(-1);
  }
  assert(mutexes.count(ml.addr));
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  // curnode().kind = ZEvent::Kind::M_DESTROY;
  // mayConflict(&ml);
  mutexes.erase(ml.addr);
  end_err();
}


void ZBuilderSC::mutex_unlock(const SymAddrSize &ml)
{
  start_err("mutex_unlock");
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
  assert(prefix_idx == (int) (prefix->size() - 1));
  if ((*prefix)[prefix_idx]->size > 1) {
    // Some invisible instructions happened before the
    // lock instruction, we split this into two events
    assert(sch_extend && !sch_replay &&
           "Lock with invisible instruction(s) should not be one event in replay trace");
    unsigned p = curnode().iid.get_pid();
    --((*prefix)[prefix_idx]->size); // Subtract the lock instruction
    // Create a new event with only the lock instruction
    ++prefix_idx;
    // Mark that thread p owns the new event
    threads[p].event_indices.push_back(prefix_idx);
    // Create the new event
    prefix->emplace_back(new ZEvent(
      IID<IPid>(IPid(p),threads[p].last_event_index()),
      threads[p].cpid,
      threads[p].executed_instructions, // first will be 1
      threads[p].executed_events, // first will be 0
      prefix->size()));
    fence(); // Each event in SC has fence
    ++threads[p].executed_events;
    assert(prefix->back()->trace_id() == (int) prefix->size() - 1);
    assert((unsigned) prefix_idx == prefix->size() - 1);
    // Graph and PO building
    graph_po_process_event(false);
  }

  assert(curnode().size == 1);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::M_LOCK;
  curnode()._observed_trace_id = mutex.last_access; // initialized with -1

  mayConflict(&ml);
  if (!sch_replay) {
    somethingToAnnotate.insert(curnode().iid.get_pid());
    threads_with_unannotated_readlock.insert(curnode().cpid());
  }
  endsWithLockFail.erase(curnode().iid.get_pid());

  assert(!mutex.locked);
  mutex.last_lock = mutex.last_access = prefix_idx;
  mutex.locked = true;
  end_err();
}


void ZBuilderSC::mutex_lock_fail(const SymAddrSize &ml) {
  start_err("mutex_lock_fail");
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
  end_err();
}


void ZBuilderSC::extendGivenTrace() {
  start_err("extendGivenTrace...");
  assert(sch_replay && !replay_trace.empty());

  std::unique_ptr<llvm::ExecutionEngine> EE(DPORDriver::create_execution_engine(M, *this, *config));

  // Run main.
  EE->runFunctionAsMain(M->getFunction("main"), {"prog"}, 0);

  // Run static destructors.
  EE->runStaticConstructorsDestructors(true);

  // Add lock event for every thread ending with a failed mutex lock attempt
  add_failed_lock_attempts();

  // Graph and PO building - process the very last event
  graph_po_process_event(true);
  assert(graph_po_constructed);

  end_err();
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
    prefix->emplace_back(new ZEvent(
      IID<IPid>(IPid(p),threads[p].last_event_index()),
      threads[p].cpid,
      threads[p].executed_instructions, // first will be 1
      threads[p].executed_events, // first will be 0
      prefix->size()));
    ++threads[p].executed_events;
    assert(prefix->back()->trace_id() == (int) prefix->size() - 1);
    assert((unsigned) prefix_idx == prefix->size() - 1);
    assert(curnode().size == 1);
    assert(curnode().kind == ZEvent::Kind::DUMMY);
    curnode().kind = ZEvent::Kind::M_LOCK;

    mayConflict(&(p_ml.second));
    // Graph and PO building
    graph_po_process_event(false);
    // Mark unannotated lock only now, this event will be checked the next time
    threads_with_unannotated_readlock.insert(curnode().cpid());
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

  for(unsigned i = 0; i < prefix->size(); ++i) {
    cmp.push_back(IID<CPid>(threads[(*prefix)[i]->iid.get_pid()].cpid,
                            (*prefix)[i]->iid.get_index()));
    cmp_md.push_back((*prefix)[i]->md);
  }
  for(unsigned i = 0; i < errors.size(); ++i){
    errs.push_back(errors[i]->clone());
  }
  Trace *t = new IIDSeqTrace(cmp,cmp_md,errs);
  t->set_blocked(false);
  return t;
}


void ZBuilderSC::graph_po_process_event(bool take_last_event)
{
  start_err("graph_po_process_event...");
  assert(prefix->size() > 1);
  assert(prefix_idx == prefix->size() - 1);
  const ZEvent * ev;
  if (take_last_event) {
    // We are at the end of interpreting, all events except
    // the very last one have already been processed
    ev = (*prefix)[prefix_idx].get();
  } else {
    // We are in the middle of interpreting, we process the second-to-last
    // event since the last one might get refuse_schedule-d
    ev = (*prefix)[prefix_idx - 1].get();
    if (graph->has_event(ev)) {
      // We are revisiting ev because some event after it got
      // refuse_schedule-d and now there is a different event after it
      // Nothing to do now since we have processed ev before
      assert(graph->events_size() == ev->trace_id() + 1);
      return;
    }
  }
  assert(graph->events_size() == ev->trace_id());
  assert(!graph->has_event(ev));

  if (prefix_idx == 1) {
    assert(!take_last_event);
    graph->inherit_lines(*po_part);
    po_full->inherit_lines(*po_part);
  }

  // Add to graph and po_full

  if (!graph->has_thread(ev->cpid())) {
    graph->add_line(ev);
    graph->add_event(ev);
    assert(!po_full->spans_thread(ev->cpid()));
    po_full->add_line(ev);
    po_full->add_event(ev);
  } else {
    graph->add_event(ev);
    assert(po_full->spans_thread(ev->cpid()));
    po_full->add_event(ev);
  }
  assert(graph->has_thread(ev->cpid()));
  assert(po_full->spans_thread(ev->cpid()));
  assert(graph->has_event(ev));
  assert(po_full->spans_event(ev));

  // Handle specific event types

  if (is_spawn(ev)) {
    spawns.push_back(ev);
    if (threads_past_annotated_region.count(ev->cpid())) {
      assert(!threads_past_annotated_region.count(ev->childs_cpid()));
      threads_past_annotated_region.insert(ev->childs_cpid());
    }
    // Cache - last spawn
    graph->_cache.last_spawn[ev->cpid()] = ev;
  }

  if (is_join(ev)) {
    joins.push_back(ev);
    if (threads_past_annotated_region.count(ev->childs_cpid()) &&
        !threads_past_annotated_region.count(ev->cpid())) {
      threads_past_annotated_region.insert(ev->cpid());
    }
  }

  if (is_write(ev)) {
    // Last write of thread
    if (!last_write_of_thread.count(ev->cpid()))
      last_write_of_thread.emplace
        (ev->cpid(), std::unordered_map<SymAddrSize, const ZEvent *>());
    last_write_of_thread[ev->cpid()][ev->ml()] = ev;
    // Cache - writes
    auto& writes = graph->_cache.writes;
    if (!writes.count(ev->ml()))
      writes.emplace
        (ev->ml(), std::unordered_map<CPid, std::vector<const ZEvent *>>());
    if (!writes[ev->ml()].count(ev->cpid()))
      writes[ev->ml()].emplace(ev->cpid(), std::vector<const ZEvent *>());
    writes[ev->ml()][ev->cpid()].push_back(ev);
  }

  if (is_read(ev)) {
    // Cache - writes
    auto& writes = graph->_cache.writes;
    if (!writes.count(ev->ml()))
      writes.emplace
        (ev->ml(), std::unordered_map<CPid, std::vector<const ZEvent *>>());
    // Cache - local_write
    auto& local_write = graph->_cache.local_write;
    assert(!local_write.count(ev));
    if (last_write_of_thread.count(ev->cpid()) &&
        last_write_of_thread[ev->cpid()].count(ev->ml())) {
      assert(is_write(last_write_of_thread.at(ev->cpid()).at(ev->ml())));
      local_write.emplace(ev, last_write_of_thread[ev->cpid()][ev->ml()]);
    }
    else
      local_write.emplace(ev, nullptr);
  }

  if (is_lock(ev)) {
    // Cache - last lock
    graph->_cache.last_lock[ev->cpid()] = ev;
  }

  if (is_unlock(ev)) {
    // Cache - last unlock
    auto& last_unlock = graph->_cache.last_unlock;
    if (!last_unlock.count(ev->ml())) {
      last_unlock.emplace(ev->ml(), std::map<CPid, const ZEvent *>());
    }
    last_unlock[ev->ml()][ev->cpid()] = ev;
    assert(last_unlock.at(ev->ml()).at(ev->cpid()) == ev);
  }

  // Add to po_part if not past annotated region

  if (!threads_past_annotated_region.count(ev->cpid())) {
    if (!po_part->spans_thread(ev->cpid())) {
      po_part->add_line(ev);
    }
    if (!po_part->spans_event(ev)) {
      assert(sch_extend || ev->is_write_of_atomic_event());
      po_part->add_event(ev);
    }
    assert(po_part->spans_event(ev));

    if (threads_with_unannotated_readlock.count(ev->cpid())) {
      // We are at the first unannotated read/lock of this thread
      assert(sch_extend);
      assert(is_read(ev) || is_lock(ev));
      threads_past_annotated_region.insert(ev->cpid());
    }
  }

  if (!take_last_event) {
    end_err("not_last");
    return;
  }

  // We processed all events

  // EDGES - spawns
  for (const ZEvent *spwn : spawns) {
    assert(graph->has_thread(spwn->childs_cpid()));
    assert(!(*graph)(spwn->childs_cpid()).empty());
    const ZEvent *nthr = graph->event(spwn->childs_cpid(), 0);

    assert(!po_full->has_edge(nthr, spwn));
    if (!po_full->has_edge(spwn, nthr))
      po_full->add_edge(spwn, nthr);

    if (!po_part->spans_event(spwn))
      continue;
    assert(po_part->spans_event(nthr));
    assert(!po_part->has_edge(nthr, spwn));
    if (!po_part->has_edge(spwn, nthr))
      po_part->add_edge(spwn, nthr);
  }

  // EDGES - joins
  for (const ZEvent *jn : joins) {
    assert(graph->has_thread(jn->childs_cpid()));
    assert(!(*graph)(jn->childs_cpid()).empty());
    int lastev_idx = (*graph)(jn->childs_cpid()).size() - 1;
    assert(lastev_idx >= 0);
    const ZEvent *wthr = graph->event(jn->childs_cpid(), lastev_idx);

    assert(!po_full->has_edge(jn, wthr));
    if (!po_full->has_edge(wthr, jn))
      po_full->add_edge(wthr, jn);

    if (!po_part->spans_event(jn))
      continue;
    assert(po_part->spans_event(wthr));
    assert(!po_part->has_edge(jn, wthr));
    if (!po_part->has_edge(wthr, jn))
      po_part->add_edge(wthr, jn);
  }

  #ifndef NDEBUG
  graph_po_constructed = true;
  #endif
  end_err("last");
}
