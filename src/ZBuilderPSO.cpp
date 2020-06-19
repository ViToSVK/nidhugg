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
#include "ZBuilderPSO.h"
#include "ZExplorer.h"


// Use at the very beginning to get an initial trace
ZBuilderPSO::ZBuilderPSO
(const Configuration &conf, llvm::Module *m, unsigned s_r_i)
  : PSOTraceBuilder(conf),
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
ZBuilderPSO::ZBuilderPSO
(const Configuration &conf, llvm::Module *m, std::vector<ZEvent>&& tr)
  : PSOTraceBuilder(conf),
    sch_replay(true), sch_extend(false),
    replay_trace(std::move(tr))
{
  config = &conf;
  M = m;
  prefix.reserve(replay_trace.size() + 64);
  assert(tr.empty() && !replay_trace.empty());
  ext_from_id = replay_trace.size();
}


bool ZBuilderPSO::reset()
{
  #ifndef NDEBUG
    std::cout << "RUNNING DEBUG VERSION\n";
  #endif

  if (this->has_error()) {
    this->error_trace = this->get_trace();
    return true;
  }

  // Lock event for every thread ending with a failed mutex lock attempt
  // Do not add in the maximum-trace exploration
  // add_failed_lock_attempts();

  // Construct the explorer with the original ZBuilder pointer
  ZExplorer explorer(*this);

  // Construct initial trace
  ZTrace initial_trace(
    nullptr, std::vector<ZEvent>(), std::vector<std::unique_ptr<ZEvent>>(),
    ZAnnotation(), std::set<ZEventID>());

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


bool ZBuilderPSO::schedule(int *proc, int *aux, int *alt, bool *dryrun)
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


bool ZBuilderPSO::schedule_replay_trace(int *proc, int *aux)
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
    if (prefix_idx == 0) {
      cpidMainToIpid.emplace(replay_trace[prefix_idx].cpid().get_proc_seq(),
                             replay_trace[prefix_idx].iid.get_pid());
    }
    const ZEvent& toReplay = replay_trace[prefix_idx];
    assert(cpidMainToIpid.count(toReplay.cpid().get_proc_seq()));
    assert(!isWriteM(toReplay) ||
           (cpidMlToIpid.count(toReplay.cpid().get_proc_seq()) &&
            cpidMlToIpid[toReplay.cpid().get_proc_seq()].count(toReplay.ml)));
    unsigned p = (isWriteM(toReplay))
      ? cpidMlToIpid[toReplay.cpid().get_proc_seq()][toReplay.ml]
      : cpidMainToIpid[toReplay.cpid().get_proc_seq()];
    // Create the new event
    assert(replay_trace[prefix_idx].instruction_id == (unsigned)
           threads[p].executed_instructions + 1 && "Inconsistent scheduling");
    assert(replay_trace[prefix_idx].event_id() ==
           threads[p].executed_events && "Inconsistent scheduling");
    // ++threads[p].executed_instructions; Do it later after this block
    prefix.emplace_back(IID<IPid>(IPid(p),prefix.size()),
                        threads[p].cpid,
                        threads[p].executed_instructions + 1, // +1 so that first will be 1
                        threads[p].executed_events, // so that first will be 0
                        prefix.size());
    // Mark that thread executes a new event
    ++threads[p].executed_events;
    assert(prefix.back().trace_id() == (int) prefix.size() - 1);
  } else {
    // Next instruction is a continuation of the current event
    assert(replay_trace[prefix_idx].size > prefix[prefix_idx].size);
    // Increase the size of the current event, it will be scheduled
    ++prefix[prefix_idx].size;
    assert(replay_trace[prefix_idx].size >=
           prefix[prefix_idx].size && "Inconsistent scheduling");
  }

  assert((unsigned) prefix_idx < replay_trace.size());
  const ZEvent& toReplay = replay_trace[prefix_idx];
  assert(cpidMainToIpid.count(toReplay.cpid().get_proc_seq()));
  assert(!isWriteM(toReplay) ||
         (cpidMlToIpid.count(toReplay.cpid().get_proc_seq()) &&
          cpidMlToIpid[toReplay.cpid().get_proc_seq()].count(toReplay.ml)));
  unsigned p = (isWriteM(toReplay))
    ? cpidMlToIpid[toReplay.cpid().get_proc_seq()][toReplay.ml]
    : cpidMainToIpid[toReplay.cpid().get_proc_seq()];
  // Mark that thread p executes a new instruction
  ++threads[p].executed_instructions;
  assert(threads[p].available);

  bool ret = schedule_thread(proc, aux, p);
  assert(ret && "Bug in scheduling: could not reproduce a given replay trace");

  return ret;
}


bool ZBuilderPSO::schedule_arbitrarily(int *proc, int *aux)
{
  assert(!sch_replay && sch_extend);

  // We prefer auxiliary before real threads

  for(int p_aux : available_auxs){ // Loop through auxiliary threads
    assert(p_aux >= 0);
    if (schedule_thread(proc, aux, (unsigned) p_aux))
      return true;
  }
  for(int p_real : available_threads){ // Loop through real threads
    assert(p_real >= 0);
    if (schedule_thread(proc, aux, (unsigned) p_real))
      return true;
  }

  // We did not schedule anything
  prefix.shrink_to_fit();
  return false;
}


// Schedule the next instruction to be the next instruction from thread p
bool ZBuilderPSO::schedule_thread(int *proc, int *aux, unsigned p)
{
  if (threads[p].available && !threads[p].sleeping &&
      (conf.max_search_depth < 0 || threads[p].clock[p] < conf.max_search_depth) &&
      (!threads[p].cpid.is_auxiliary() || is_aux_at_head(p))) {

    if (sch_extend) {
      update_prefix(p);
    } else {
      assert(sch_replay);
    }

    // set the scheduled thread and aux
    *proc = threads[p].proc;
    *aux = (threads[p].cpid.is_auxiliary())
      ? threads[p].cpid.get_aux_index() : -1;

    return true;
  }

  return false;
}


// Used only when we are not replaying replay_trace
void ZBuilderPSO::update_prefix(unsigned p)
{
  assert(!sch_replay && sch_extend);
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
    // Create the new event
    prefix.emplace_back(IID<IPid>(IPid(p),prefix.size()),
                        threads[p].cpid,
                        threads[p].executed_instructions, // first will be 1
                        threads[p].executed_events, // first will be 0
                        prefix.size());
    ++threads[p].executed_events;
    assert(prefix.back().trace_id() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }
}


/* ******************************************************* */
/* ADDRESSING FEEDBACK FROM INTERPRETER AFTER SCHEDULING   */
/* ******************************************************* */


void ZBuilderPSO::refuse_schedule()
{
  //llvm::errs() << " REFUSESCH";
  assert(prefix_idx == int(prefix.size())-1);
  assert(!prefix.back().may_conflict);

  unsigned p = curnode().iid.get_pid();
  assert(!threads[p].cpid.is_auxiliary() && "Only real threads can be refused");

  // Here used to be:
  // --threads[p].event_indices[p];

  --threads[p].executed_instructions; //

  if (curnode().size == 1) {
    // Refused instruction wanted to create a new event,
    // therefore this entire event needs to be deleted
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
  mark_unavailable_ipid(p);
}


void ZBuilderPSO::mayConflict(const SymAddrSize *ml)
{
  auto& curn = curnode();
  curn.may_conflict = true;
  if (ml) curn.ml = *ml;

  #ifndef NDEBUG
  bool consistent = (!sch_replay ||
                     replay_trace[prefix_idx].kind == prefix[prefix_idx].kind);
  if (!consistent) {
    llvm::errs() << "TRACE_TO_REPLAY\n";
    dump_trace(replay_trace);
    llvm::errs() << "\nACTUALLY_REPLAYED\n";
    dump_trace(prefix);
    llvm::errs() << "Problematic event (didn't go any further):\n";
    llvm::errs() << "TRACE_TO_REPLAY[" << prefix_idx << "]   ::: ";
    replay_trace[prefix_idx].dump();
    llvm::errs() << "ACTUALLY_REPLAYED[" << prefix_idx << "] ::: ";
    prefix[prefix_idx].dump();
  }
  assert(consistent);
  #endif
}


void ZBuilderPSO::metadata(const llvm::MDNode *md)
{
  if (md) curnode().md = md;
}


IID<CPid> ZBuilderPSO::get_iid() const
{
  IPid pid = curnode().iid.get_pid();
  int idx = curnode().iid.get_index();
  return IID<CPid>(threads[pid].cpid,idx);
}


void ZBuilderPSO::spawn()
{
  //llvm::errs() << " SPAWN";
  assert(!dryrun);
  IPid parent_ipid = curnode().iid.get_pid();
  assert(!threads[parent_ipid].cpid.is_auxiliary());
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::SPAWN;
  mayConflict();

  IPid child_ipid = threads.size();
  CPid child_cpid = CPS.spawn(threads[parent_ipid].cpid);
  assert(!child_cpid.is_auxiliary());
  cpidMainToIpid.emplace(child_cpid.get_proc_seq(), child_ipid);
  curnode().childs_cpid = child_cpid;
  int proc = 0;
  for(unsigned i = 0; i < threads.size(); ++i){
    proc = std::max(proc,threads[i].proc+1);
  }

  proc_to_ipid.push_back(child_ipid);
  threads.push_back(Thread(proc,child_cpid,threads[parent_ipid].clock,parent_ipid));
  mark_available_ipid(child_ipid);
  spawned_something = true;
}


void ZBuilderPSO::join(int tgt_proc)
{
  //llvm::errs() << " JOIN";
  assert(!dryrun);
  assert(0 <= tgt_proc && tgt_proc < int(proc_to_ipid.size()));
  unsigned p = curnode().iid.get_pid();
  assert(!threads[p].cpid.is_auxiliary());
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  // We make sure that every join is an event with exactly
  // one instruction - only the join instruction itself
  assert(prefix_idx == (int) (prefix.size() - 1));
  if (prefix[prefix_idx].size > 1) {
    // Some invisible instructions happened before the
    // join instruction, we split this into two events
    --prefix[prefix_idx].size; // Subtract the join instruction
    // Create a new event with only the join instruction
    ++prefix_idx;
    // Create the new event
    prefix.emplace_back(IID<IPid>(IPid(p),prefix.size()),
                        threads[p].cpid,
                        threads[p].executed_instructions, // first will be 1
                        threads[p].executed_events, // first will be 0
                        prefix.size());
    ++threads[p].executed_events;
    assert(prefix.back().trace_id() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }

  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::JOIN;
  IPid tgt_ipid = proc_to_ipid[tgt_proc];
  curnode().childs_cpid = threads[tgt_ipid].cpid;
  mayConflict();
  //curnode().clock += threads[tgt_ipid].clock;
  //curnode().clock += threads[tgt_ipid].aux_clock_sum;
  //threads[ipid].clock += curnode().clock;
}


void ZBuilderPSO::store(const SymData &sd, int val)
{
  assert(!dryrun);
  unsigned ipid = curnode().iid.get_pid();
  assert(!threads[ipid].cpid.is_auxiliary());

  const SymAddrSize& ml = sd.get_ref();
  // Stores to local memory on stack may not conflict
  assert(ml.addr.block.is_global() || ml.addr.block.is_heap());

  // Put the store to buffer
  for(SymAddr b : ml){
    threads[ipid].store_buffers[b].push_back(PendingStoreByte(ml,threads[ipid].clock,last_md));
  }
  IPid upd_ipid;
  auto it = threads[ipid].byte_to_aux.find(ml.addr);
  if(it == threads[ipid].byte_to_aux.end()){
    /* Create new auxiliary thread */
    int aux_idx = int(threads[ipid].aux_to_byte.size());
    upd_ipid = int(threads.size());
    CPid upd_cpid = CPS.new_aux(threads[ipid].cpid);
    assert(upd_cpid.is_auxiliary());
    threads.push_back(Thread(threads[ipid].proc,upd_cpid,threads[ipid].clock,ipid));
    threads[ipid].byte_to_aux[ml.addr] = aux_idx;
    threads[ipid].aux_to_byte.push_back(ml.addr);
    threads[ipid].aux_to_ipid.push_back(upd_ipid);
    if (!cpidMlToIpid.count(upd_cpid.get_proc_seq()))
      cpidMlToIpid.emplace(upd_cpid.get_proc_seq(),
                           std::unordered_map<SymAddrSize, unsigned>());
    assert(!cpidMlToIpid[upd_cpid.get_proc_seq()].count(ml));
    cpidMlToIpid[upd_cpid.get_proc_seq()].emplace(ml, upd_ipid);
  }else{
    upd_ipid = threads[ipid].aux_to_ipid[it->second];
  }
  mark_available_ipid(upd_ipid);

  // Global stores happening with just one thread existing may not conflict
  // but we have to record them as such anyway to keep track of them
  //llvm::errs() << " BUFFERWRITE_" << ml.to_string() << " ";
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::WRITEB;
  curnode().value = val;
  mayConflict(&ml);

  // Enqueue the update into visible store queue
  if (!visibleStoreQueue.count(ipid))
    visibleStoreQueue.emplace(ipid, std::unordered_map<SymAddrSize, std::list<int>>());
  if (!visibleStoreQueue[ipid].count(ml))
    visibleStoreQueue[ipid].emplace(ml, std::list<int>());
  visibleStoreQueue[ipid][ml].push_back(prefix_idx);
}


void ZBuilderPSO::atomic_store(const SymData &sd)
{
  assert(!dryrun);

  IPid uipid = curnode().iid.get_pid(); // ID of the thread changing the memory
  assert(threads[uipid].cpid.is_auxiliary() &&
         "Only auxiliary threads can perform memory-writes");
  IPid tipid = threads[uipid].parent; // ID of the (real) thread that issued the store

  const SymAddrSize& ml = sd.get_ref();
  // We consider all Stack writes as invisible and atomic,
  // hence we do not put them into the store buffer
  assert((ml.addr.block.is_global() || ml.addr.block.is_heap()) &&
         "No Stack writes are allowed in the store buffer");

  assert(threads[tipid].store_buffers.count(ml.addr));
  assert(threads[tipid].store_buffers[ml.addr].size());
  assert(threads[tipid].store_buffers[ml.addr].front().ml == ml);
  const PendingStoreByte &pst = threads[tipid].store_buffers[ml.addr].front();
  //curnode().clock += pst.clock;
  //threads[uipid].clock += pst.clock;
  //curnode().origin_iid = IID<IPid>(tipid,pst.clock[tipid]);
  curnode().md = pst.md;

  /* Register in memory */
  int last_rowe = threads[tipid].store_buffers[ml.addr].front().last_rowe;
  for(SymAddr b : ml){
    ByteInfo &bi = mem[b];
    bi.last_update = prefix_idx;
    bi.last_update_ml = ml;
    if(0 <= last_rowe){
      bi.last_read[threads[tipid].proc] = last_rowe;
    }
  }

  // Remove pending store from buffer
  for(SymAddr b : ml){
    std::vector<PendingStoreByte> &sb = threads[tipid].store_buffers[b];
    for(unsigned i = 0; i < sb.size() - 1; ++i){
      sb[i] = sb[i+1];
    }
    sb.pop_back();
    if(sb.empty()){
      threads[tipid].store_buffers.erase(b);
      mark_unavailable_ipid(uipid);
    }
  }
  //threads[tipid].aux_clock_sum += curnode().clock;

  //llvm::errs() << " MEMORYWRITE_" << ml.to_string() << " ";
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::WRITEM;
  mayConflict(&ml);

  // Locate corresponding buffer-write
  assert(visibleStoreQueue.count(tipid) &&
         visibleStoreQueue[tipid].count(ml) &&
         !visibleStoreQueue[tipid][ml].empty());
  curnode().write_other_trace_id = visibleStoreQueue[tipid][ml].front();
  assert(curnode().write_other_trace_id >= 0 &&
         curnode().write_other_trace_id <= (int) prefix.size() - 2);
  assert(prefix[curnode().write_other_trace_id].kind == ZEvent::Kind::WRITEB);
  assert(prefix[curnode().write_other_trace_id].ml == ml &&
         "Inconsistent memory locations");
  curnode().value = prefix[curnode().write_other_trace_id].value;
  assert(prefix[curnode().write_other_trace_id].write_other_trace_id == -1);
  prefix[curnode().write_other_trace_id].write_other_trace_id = prefix_idx;

  // Dequeue oldest update from queue
  visibleStoreQueue[tipid][ml].pop_front();

  lastWrite[ml] = prefix_idx;
}


void ZBuilderPSO::load(const SymAddrSize &ml, int val)
{
  // Loads from local memory on stack may not conflict
  // Also global loads happening with just one thread existing may not conflict
  assert(!dryrun);
  assert(!threads[curnode().iid.get_pid()].cpid.is_auxiliary());
  if ((ml.addr.block.is_global() || ml.addr.block.is_heap()) &&
      spawned_something) {
    //llvm::errs() << " READ_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << " ";
    assert(curnode().kind == ZEvent::Kind::DUMMY);
    curnode().kind = ZEvent::Kind::READ;
    curnode().value = val;
    mayConflict(&ml);

    // Set ID of observed event
    // last memory-write
    int obs_idx = (lastWrite.count(ml)) ? lastWrite[ml] : -1;
    // check for buffer-writes present in queue
    unsigned p = curnode().iid.get_pid();
    if (!visibleStoreQueue.count(p))
      visibleStoreQueue.emplace
        (p, std::unordered_map<SymAddrSize, std::list<int>>());
    if (!visibleStoreQueue[p].count(ml))
      visibleStoreQueue[p].emplace(ml, std::list<int>());
    else if (!visibleStoreQueue[p][ml].empty())
      obs_idx = visibleStoreQueue[p][ml].back();
    curnode().observed_trace_id = obs_idx;
  }
}


void ZBuilderPSO::fence()
{
  assert(!dryrun);
  #ifndef NDEBUG
  IPid ipid = curnode().iid.get_pid();
  assert(!threads[ipid].cpid.is_auxiliary());
  assert(threads[ipid].all_buffers_empty());
  if (visibleStoreQueue.count(ipid) &&
      !visibleStoreQueue[ipid].empty()) {
    for (const auto& ml_list : visibleStoreQueue[ipid])
      assert(visibleStoreQueue[ipid][ml_list.first].empty());
  }
  #endif
  curnode().fence = true;
}


void ZBuilderPSO::mutex_init(const SymAddrSize &ml)
{
  //llvm::errs() << " M_INIT_" << ml.to_string() << " ";
  assert(!dryrun);
  assert(!threads[curnode().iid.get_pid()].cpid.is_auxiliary());
  fence();
  assert(!mutexes.count(ml.addr));
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::M_INIT;
  mayConflict(&ml);
  mutexes[ml.addr] = Mutex(-1); // prefix_idx
}


void ZBuilderPSO::mutex_destroy(const SymAddrSize &ml)
{
  //llvm::errs() << " M_DESTROY_" << ml.to_string() << " ";
  assert(!dryrun);
  assert(!threads[curnode().iid.get_pid()].cpid.is_auxiliary());
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


void ZBuilderPSO::mutex_unlock(const SymAddrSize &ml)
{
  //llvm::errs() << " M_UNLOCK_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << " ";
  assert(!dryrun);
  assert(!threads[curnode().iid.get_pid()].cpid.is_auxiliary());
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


void ZBuilderPSO::mutex_lock(const SymAddrSize &ml)
{
  //llvm::errs() << " M_LOCK_" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << "_";
  assert(!dryrun);
  int ipid = curnode().iid.get_pid();
  assert(!threads[ipid].cpid.is_auxiliary());
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
    // Create the new event
    prefix.emplace_back(IID<IPid>(IPid(p),prefix.size()),
                        threads[p].cpid,
                        threads[p].executed_instructions, // first will be 1
                        threads[p].executed_events, // first will be 0
                        prefix.size());
    ++threads[p].executed_events;
    assert(prefix.back().trace_id() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
  }

  assert(curnode().size == 1);
  assert(curnode().kind == ZEvent::Kind::DUMMY);
  curnode().kind = ZEvent::Kind::M_LOCK;
  curnode().observed_trace_id = mutex.last_access; // initialized with -1

  mayConflict(&ml);
  endsWithLockFail.erase(ipid);

  assert(!mutex.locked);
  mutex.last_lock = mutex.last_access = prefix_idx;
  mutex.locked = true;
}


void ZBuilderPSO::mutex_lock_fail(const SymAddrSize &ml) {
  //llvm::errs() << " M_LOCKFAIL" << ml.to_string() << "_ipid:" << curnode().iid.get_pid() << "_";
  assert(!dryrun);
  int ipid = curnode().iid.get_pid();
  assert(!threads[ipid].cpid.is_auxiliary());
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

  endsWithLockFail.emplace(ipid, ml);
}


ZTraceExtension ZBuilderPSO::extendGivenTrace() {
  assert(sch_replay && !replay_trace.empty());

  std::unique_ptr<llvm::ExecutionEngine> EE(DPORDriver::create_execution_engine(M, *this, *config));

  // Run main.
  EE->runFunctionAsMain(M->getFunction("main"), {"prog"}, 0);

  // Run static destructors.
  EE->runStaticConstructorsDestructors(true);

  // Lock event for every thread ending with a failed mutex lock attempt
  // Do not add in the maximum-trace exploration
  // add_failed_lock_attempts();

  ZTraceExtension res(
    std::move(prefix), ext_from_id,
    someThreadAssumeBlocked, !endsWithLockFail.empty());
  assert(prefix.empty());
  return res;
}


void ZBuilderPSO::add_failed_lock_attempts() {
  // Add lock event for every thread ending with a failed mutex lock attempt
  for (auto p_ml : endsWithLockFail) {
    unsigned p = p_ml.first;
    assert(!threads[p].cpid.is_auxiliary());
    ++threads[p].executed_instructions;
    ++prefix_idx;
    // Create the new event
    prefix.emplace_back(IID<IPid>(IPid(p),prefix.size()),
                        threads[p].cpid,
                        threads[p].executed_instructions, // first will be 1
                        threads[p].executed_events, // first will be 0
                        prefix.size());
    ++threads[p].executed_events;
    assert(prefix.back().trace_id() == (int) prefix.size() - 1);
    assert((unsigned) prefix_idx == prefix.size() - 1);
    assert(curnode().size == 1);
    assert(curnode().kind == ZEvent::Kind::DUMMY);
    curnode().kind = ZEvent::Kind::M_LOCK;

    fence();
    mayConflict(&(p_ml.second));
  }
}


Trace *ZBuilderPSO::get_trace() const
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
