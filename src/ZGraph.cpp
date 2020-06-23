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

#include "ZHelpers.h"
#include "ZGraph.h"
#include "ZBuilderTSO.h"

static const bool DEBUG = false;
#include "ZDebug.h"


// Empty
ZGraph::ZGraph(MemoryModel model)
  : model(model),
    init(ZEvent(true)),
    lines(),
    proc_seq_to_spawn(),
    thread_aux_to_line_id(),
    threads_auxes(),
    proc_seq_to_thread_id(),
    thr_to_lin_id(),
    lin_to_thr_id(),
    event_to_position(),
    po(*this),
    cache()
{
  assert(empty());
}


// Moving
ZGraph::ZGraph(ZGraph&& oth)
  : model(oth.model),
    init(std::move(oth.init)),
    lines(std::move(oth.lines)),
    proc_seq_to_spawn(std::move(oth.proc_seq_to_spawn)),
    thread_aux_to_line_id(std::move(oth.thread_aux_to_line_id)),
    threads_auxes(std::move(oth.threads_auxes)),
    proc_seq_to_thread_id(std::move(oth.proc_seq_to_thread_id)),
    thr_to_lin_id(std::move(oth.thr_to_lin_id)),
    lin_to_thr_id(std::move(oth.lin_to_thr_id)),
    event_to_position(std::move(oth.event_to_position)),
    po(std::move(oth.po)),
    cache(std::move(oth.cache))
{
  assert(oth.empty());
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZGraph& gr)
{
  out << gr.getPo().to_string();
  return out;
}


/* *************************** */
/* LINES                       */
/* *************************** */

const LineT& ZGraph::operator()(std::pair<unsigned, int> ids) const
{
  assert(ids.first != INT_MAX && "Called for initial event");
  auto it = thread_aux_to_line_id.find(ids);
  assert(it != thread_aux_to_line_id.end());
  return lines[it->second];
}


const LineT& ZGraph::operator()(unsigned thread_id, int aux_id) const
{
  return operator()(std::pair<unsigned, int>(thread_id, aux_id));
}


const ZEvent * const ZGraph::getEvent(unsigned thread_id, int aux_id, unsigned event_id) const
{
  assert(thread_id != INT_MAX && "Called for initial event");
  assert(hasThreadAux(thread_id, aux_id));
  const LineT& line = this->operator()(thread_id, aux_id);
  assert(event_id < line.size());
  assert(hasEvent(line[event_id]));
  return line[event_id];
}


const ZEvent * const ZGraph::getEvent(const ZEventID& id) const
{
  assert(id != ZEventID(true) && "Called for initial event");
  unsigned thr = proc_seq_to_thread_id.at(id.cpid().get_proc_seq());
  int aux = id.cpid().is_auxiliary() ? id.cpid().get_aux_index() : -1;
  unsigned ev = id.event_id();
  return getEvent(thr, aux, ev);
}


void ZGraph::addLine(const ZEvent *ev)
{
  assert(ev);
  assert(!hasEvent(ev));
  assert(ev->thread_id() < 20 && "Thread ID not set up yet");
  assert(!hasThreadAux(ev->thread_id(), ev->aux_id()));

  auto key = std::pair<unsigned, int>(ev->thread_id(), ev->aux_id());
  thread_aux_to_line_id.emplace(key, lines.size());
  if (!threads_auxes.count(ev->thread_id()))
    threads_auxes.emplace(ev->thread_id(), std::set<int>());
  assert(!threads_auxes.at(ev->thread_id()).count(ev->aux_id()));
  threads_auxes[ev->thread_id()].emplace(ev->aux_id());
  lines.push_back(LineT());
  lines.back().reserve(16);
}


void ZGraph::addEvent(const ZEvent *ev)
{
  assert(ev);
  assert(!hasEvent(ev));
  assert(ev->thread_id() < 20 && "Thread ID not set up yet");
  assert(hasThreadAux(ev->thread_id(), ev->aux_id()));

  auto key = std::pair<unsigned, int>(ev->thread_id(), ev->aux_id());
  auto it = thread_aux_to_line_id.find(key);
  assert(it != thread_aux_to_line_id.end());
  LineT& line = lines[it->second];

  unsigned event_id = line.size();
  assert(event_id == ev->event_id());
  assert(model == MemoryModel::SC || model == MemoryModel::TSO ||
         ev->aux_id() == -1 || line.empty() || sameMl(line[0], ev));
  event_to_position.emplace(ev, std::pair<unsigned,unsigned>(it->second, event_id));
  line.push_back(ev);
  assert(hasEvent(ev));
}


void ZGraph::shrink()
{
  lines.shrink_to_fit();
  for (unsigned i=0; i<lines.size(); ++i)
    lines[i].shrink_to_fit();
  for (auto& ml_x : cache.wm)
    for (auto& thr_mw : ml_x.second)
      thr_mw.second.shrink_to_fit();
}


unsigned ZGraph::lineID(unsigned thread_id, int aux_id) const
{
  assert(thread_id != INT_MAX && "Called for initial event");
  assert(hasThreadAux(thread_id, aux_id));
  auto key = std::pair<unsigned,int>(thread_id, aux_id);
  auto it = thread_aux_to_line_id.find(key);
  assert(it != thread_aux_to_line_id.end());
  return it->second;
}


unsigned ZGraph::lineID(const ZEvent *ev) const
{
  assert(ev);
  assert(!isInitial(ev) && "Called for initial event");
  return lineID(ev->thread_id(), ev->aux_id());
}


bool ZGraph::hasThreadAux(std::pair<unsigned, int> ids) const
{
  assert(ids.first != INT_MAX && "Called for initial event");
  assert(!threads_auxes.count(ids.first) ||
         !threads_auxes.at(ids.first).count(ids.second) ||
         thread_aux_to_line_id.count(ids));
  assert(!thread_aux_to_line_id.count(ids) ||
         threads_auxes.at(ids.first).count(ids.second));
  return thread_aux_to_line_id.count(ids);
}


bool ZGraph::hasThreadAux(unsigned thread_id, int aux_id) const
{
  return hasThreadAux(std::pair<unsigned, int>(thread_id, aux_id));
}


unsigned ZGraph::number_of_threads() const
{
  return threads_auxes.size();
}


std::set<unsigned> ZGraph::get_threads() const
{
  std::set<unsigned> res;
  for (const auto& tid_auxes : threads_auxes) {
    assert(tid_auxes.first >= 0);
    res.emplace((unsigned) tid_auxes.first);
  }
  return res;
}


const std::set<int>& ZGraph::auxes(unsigned thread_id) const
{
  assert(thread_id != INT_MAX && "Called for initial event");
  assert(threads_auxes.count(thread_id));
  return threads_auxes.at(thread_id);
}


int ZGraph::auxForMl(const SymAddrSize& ml, unsigned thr) const
{
  assert(threads_auxes.count(thr));
  auto axs = auxes(thr);
  assert(!axs.empty());
  if (axs.size() == 1) {
    // No write events in this thread
    return -1;
  }
  if (axs.size() == 2 && model != MemoryModel::PSO) {
    #ifndef NDEBUG
    auto it = axs.begin();
    while (*it == -1) {
      ++it;
      assert(it != axs.end());
    }
    assert(*it == 0);
    #endif
    return 0;
  }
  // PSO below
  assert(model == MemoryModel::PSO);
  for (auto aux : axs) {
    if (aux != -1) {
      assert(!((*this)(thr, aux).empty()));
      const ZEvent *first = (*this)(thr, aux)[0];
      assert(isWriteM(first));
      if (first->ml == ml) {
        assert(first->aux_id() == aux);
        return aux;
      }
    }
  }
  // This thread has no event for this ml
  return -1;
}


void ZGraph::addThreadID(const ZEvent *ev)
{
  assert(ev);
  assert(!isInitial(ev) && "Called for initial event");
  auto it = proc_seq_to_thread_id.find(ev->cpid().get_proc_seq());
  if (it == proc_seq_to_thread_id.end()) {
    proc_seq_to_thread_id.emplace_hint(
      it, ev->cpid().get_proc_seq(), ev->thread_id());
    assert(ev->thread_id() >= 0);
  }
}


// <threadID, added_with_this_call?>
std::pair<unsigned, bool> ZGraph::getThreadID(const ZEvent *ev)
{
  assert(ev);
  assert(!isInitial(ev) && "Called for initial event");
  auto it = proc_seq_to_thread_id.find(ev->cpid().get_proc_seq());
  if (it == proc_seq_to_thread_id.end()) {
    proc_seq_to_thread_id.emplace_hint(
      it, ev->cpid().get_proc_seq(), ev->thread_id());
    assert(ev->thread_id() >= 0);
    return {(unsigned) ev->thread_id(), true};
  }
  return {it->second, false};
}


unsigned ZGraph::getThreadIDnoAdd(const std::vector<int>& proc_seq) const
{
  auto it = proc_seq_to_thread_id.find(proc_seq);
  if (it == proc_seq_to_thread_id.end())
    return 1337;
  else
    return it->second;
}


unsigned ZGraph::getThreadIDnoAdd(const ZEvent * ev) const
{
  assert(ev);
  assert(!isInitial(ev) && "Called for initial event");
  return getThreadIDnoAdd(ev->cpid().get_proc_seq());
}


std::unordered_set<std::vector<int>> ZGraph::all_proc_seq() const
{
  std::unordered_set<std::vector<int>> res;
  for (const auto& proc_idx : proc_seq_to_thread_id) {
    assert(!res.count(proc_idx.first));
    res.emplace(proc_idx.first);
  }
  return res;
}


bool ZGraph::hasEvent(const ZEvent *ev) const
{
  assert(ev);
  if (ev == initial())
    return true;
  auto it = event_to_position.find(ev);
  if (it == event_to_position.end())
    return false;
  assert(lines[it->second.first][it->second.second] == ev);
  assert(lineID(ev) == it->second.first);
  assert(ev->event_id() == it->second.second);
  return true;
}


bool ZGraph::hasEvent(const ZEventID& id) const
{
  if (id == ZEventID(true))
    return true;
  if (!proc_seq_to_thread_id.count(id.cpid().get_proc_seq()))
    return false;
  unsigned thr = proc_seq_to_thread_id.at(id.cpid().get_proc_seq());
  int aux = id.cpid().is_auxiliary() ? id.cpid().get_aux_index() : -1;
  assert(thr != INT_MAX && "Called for initial event");
  if (!hasThreadAux(thr, aux))
    return false;
  const LineT& line = this->operator()(thr, aux);
  if (id.event_id() >= line.size())
    return false;
  assert(line[id.event_id()]->id() == id);
  return true;
}


/* *************************** */
/* GRAPH CONSTRUCTION          */
/* *************************** */

void ZGraph::construct
(const std::vector<std::unique_ptr<ZEvent>>& trace,
 int pre_tau_limit, std::set<int> causes_after)
{
  start_err("starting trace-to-po");
  assert(!trace.empty());
  assert(!constructed);
  constructed = true;

  std::unordered_map<std::pair<unsigned, int>, unsigned>
    cur_evidx;
  cur_evidx.reserve(8);

  std::unordered_map
    <unsigned, std::unordered_map
     <SymAddrSize, std::list<const ZEvent *>>> store_buffer;

  std::unordered_map
    <unsigned, std::unordered_map
     <SymAddrSize, const ZEvent *>> last_mwrite;

  std::unordered_map
    <unsigned, std::unordered_map
     <SymAddrSize, const ZEvent *>> last_bwrite;

  std::vector<const ZEvent *> spawns;
  spawns.reserve(8);

  std::vector<const ZEvent *> joins;
  joins.reserve(8);

  std::unordered_map<SymAddrSize, const ZEvent *> mutex_inits;
  mutex_inits.reserve(8);

  std::unordered_map<SymAddrSize, const ZEvent *> mutex_destroys;
  mutex_destroys.reserve(8);

  std::unordered_map
    <SymAddrSize, std::unordered_map<unsigned, const ZEvent *>> mutex_first;
  mutex_first.reserve(8);

  std::unordered_map
    <SymAddrSize, std::unordered_map<unsigned, const ZEvent *>> mutex_last;
  mutex_first.reserve(8);

  for (int i = 0; i < trace.size(); ++i) {
    // Limit to events until including pre_tau, and causes_after
    if (i > pre_tau_limit && !causes_after.count(i))
      continue;
    const ZEvent * const ev = trace[i].get();
    assert(i != pre_tau_limit || isRead(ev) || isLock(ev));

    /*
    unsigned thr_idx = INT_MAX;

    // Check if this process is already known
    auto ipidit = ipid_to_thraux.find(ev->iid.get_pid());
    if (model != MemoryModel::PSO && ipidit != ipid_to_thraux.end()) {
      thr_idx = ipidit->second.first;
      assert(ev->aux_id() == ipidit->second.second);
    } else {
      thr_idx = getThreadID(ev).first;
      // add to ipid cache for faster lookup next time
      ipid_to_thraux.emplace(ev->iid.get_pid(),
                             std::pair<unsigned, int>(thr_idx,
                                                      ev->aux_id()));
    }

    ev->_thread_id = thr_idx;
    */
    addThreadID(ev);

    // Set linearization_thr_id
    assert(ev->thread_id() < 100);
    assert(thr_to_lin_id.size() == lin_to_thr_id.size());
    if (!thr_to_lin_id.count(ev->thread_id())) {
      unsigned lin_id = lin_to_thr_id.size();
      assert(!lin_to_thr_id.count(lin_id));
      thr_to_lin_id.emplace(ev->thread_id(), lin_id);
      lin_to_thr_id.emplace(lin_id, ev->thread_id());
    }

    // THIS MIGHT NOT BE NEEDED FOR MAXIMAL TRACES
    /*
    if (model == MemoryModel::PSO && isWriteM(ev)) {
      // Handle aux for pso
      int aux_idx = psoGetAux(ev);
      ev->_aux_id = aux_idx;
    }
    */

    auto thraux = std::pair<unsigned, int>(ev->thread_id(), ev->aux_id());

    // Memory-writes cannot be added if their buffer-write was not added
    if (isWriteM(ev)) {
      assert(ev->aux_id() != -1 && "Only auxiliary threads");
      assert(hasEvent(ev->write_other_ptr));
    }

    // We will add the event
    // First, add line if not already present
    if (!hasThreadAux(thraux)) {
      // A new thread/aux appeared
      assert(!cur_evidx.count(thraux));
      cur_evidx.emplace(thraux, 0);
      assert(ev->event_id() == 0);

      addLine(ev);
      po.addLine(ev);
    }
    if (!cur_evidx.count(thraux)) {
      assert(ev->event_id() == 0);
      cur_evidx.emplace(thraux, 0);
    }

    unsigned ev_idx = cur_evidx[thraux];
    assert(ev->event_id() == ev_idx);
    assert(ev_idx == (*this)(thraux).size());
    // New event
    addEvent(ev);
    po.addEvent(ev);
    ++cur_evidx[thraux];

    // Handle Fence
    if (ev->fence) {
      if (store_buffer.count(ev->thread_id()))
        for (const auto& ml_last : last_mwrite[ev->thread_id()]) {
          auto last = ml_last.second;
          assert(isWriteM(last));
          assert(!po.hasEdge(ev, last));
          if (!po.hasEdge(last, ev))
            po.addEdge(last, ev);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (isLock(ev) || isUnlock(ev)) {
      llvm::errs() << "No support for locks atm";
      abort();
    }

    // Handle Spawn
    if (isSpawn(ev)) {
      spawns.push_back(ev);
      assert(!proc_seq_to_spawn.count(ev->childs_cpid.get_proc_seq()));
      proc_seq_to_spawn.emplace(ev->childs_cpid.get_proc_seq(), ev);
      assert(number_of_threads() > 0);
      assert(ev->fence);
    }
    // Handle Join
    if (isJoin(ev))
      joins.push_back(ev);

    // Handle Buffer-Write
    if (isWriteB(ev)) {
      assert(ev->aux_id() == -1 && "Only real threads");
      // Store buffer
      if (!store_buffer.count(ev->thread_id()))
        store_buffer.emplace
          (ev->thread_id(),
           std::unordered_map
           <SymAddrSize, std::list<const ZEvent *>>());
      if (!store_buffer[ev->thread_id()].count(ev->ml))
        store_buffer[ev->thread_id()].emplace(ev->ml,
                                             std::list<const ZEvent *>());
      store_buffer[ev->thread_id()][ev->ml].push_back(ev);
      // Last Bwrite
      if (!last_bwrite.count(ev->thread_id()))
        last_bwrite.emplace
          (ev->thread_id(), std::unordered_map<SymAddrSize, const ZEvent *>());
      auto it = last_bwrite[ev->thread_id()].find(ev->ml);
      if (it != last_bwrite[ev->thread_id()].end())
        it = last_bwrite[ev->thread_id()].erase(it);
      last_bwrite[ev->thread_id()].emplace_hint
        (it, ev->ml, ev);
    }

    // Handle Memory-Write
    if (isWriteM(ev)) {
      assert(ev->aux_id() != -1 && "Only auxiliary threads");
      // Store buffer
      assert(store_buffer.count(ev->thread_id()) &&
             store_buffer[ev->thread_id()].count(ev->ml));
      // Done already during extend
      /*
      assert(!ev->write_other_ptr);
      ev->write_other_ptr = store_buffer[ev->thread_id()][ev->ml].front();
      assert(isWriteB(ev->write_other_ptr) && sameMl(ev, ev->write_other_ptr));
      assert(!ev->write_other_ptr->write_other_ptr);
      ev->write_other_ptr->write_other_ptr = ev;
      assert(ev->write_other_ptr->trace_id() == ev->write_other_trace_id);
      assert(ev->write_other_ptr->write_other_trace_id == ev->trace_id());
      */
      store_buffer[ev->thread_id()][ev->ml].pop_front();
      // WB -> WM thread order
      assert(!po.hasEdge(ev, ev->write_other_ptr));
      if (!po.hasEdge(ev->write_other_ptr, ev))
        po.addEdge(ev->write_other_ptr, ev);
      // Last Mwrite
      if (!last_mwrite.count(ev->thread_id()))
        last_mwrite.emplace
          (ev->thread_id(), std::unordered_map<SymAddrSize, const ZEvent *>());
      auto it = last_mwrite[ev->thread_id()].find(ev->ml);
      if (it != last_mwrite[ev->thread_id()].end())
        it = last_mwrite[ev->thread_id()].erase(it);
      last_mwrite[ev->thread_id()].emplace_hint
        (it, ev->ml, ev);
      // Cache - wm
      if (!cache.wm.count(ev->ml))
        cache.wm.emplace
          (ev->ml, std::unordered_map<unsigned, std::vector<const ZEvent *>>());
      if (!cache.wm[ev->ml].count(ev->thread_id()))
        cache.wm[ev->ml].emplace(ev->thread_id(), std::vector<const ZEvent *>());
      cache.wm[ev->ml][ev->thread_id()].push_back(ev);
    }

    // Handle Read
    if (isRead(ev)) {
      // Cache - wm
      if (!cache.wm.count(ev->ml))
        cache.wm.emplace
          (ev->ml, std::unordered_map<unsigned, std::vector<const ZEvent *>>());
      // Cache - readWB
      assert(!cache.readWB.count(ev));
      if (last_bwrite.count(ev->thread_id()) &&
          last_bwrite[ev->thread_id()].count(ev->ml)) {
        assert(isWriteB(last_bwrite[ev->thread_id()][ev->ml]));
        cache.readWB.emplace(ev, last_bwrite[ev->thread_id()][ev->ml]);
      }
      else
        cache.readWB.emplace(ev, nullptr);
    }

    // Handle Mutex events
    if (isMutexInit(ev)) {
      assert(!mutex_inits.count(ev->ml));
      mutex_inits.emplace(ev->ml, ev);
    }
    if (isMutexDestroy(ev)) {
      assert(!mutex_destroys.count(ev->ml));
      mutex_destroys.emplace(ev->ml, ev);
    }
    if (isLock(ev) || isUnlock(ev)) {
      // For each ml, for each thread, remember first access
      if (!mutex_first.count(ev->ml)) {
        mutex_first.emplace(ev->ml, std::unordered_map<unsigned, const ZEvent *>());
        mutex_first[ev->ml].reserve(8);
      }
      if (!mutex_first[ev->ml].count(ev->thread_id()))
        mutex_first[ev->ml].emplace(ev->thread_id(), ev);
      // For each ml, for each thread, remember last access
      if (!mutex_last.count(ev->ml)) {
        mutex_last.emplace(ev->ml, std::unordered_map<unsigned, const ZEvent *>());
        mutex_last[ev->ml].reserve(8);
      }
      mutex_last[ev->ml][ev->thread_id()] = ev;
    }
    if (isUnlock(ev)) {
      // Cache - unl
      if (!cache.unl.count(ev->ml))
        cache.unl.emplace
          (ev->ml, std::unordered_map<unsigned, std::vector<const ZEvent *>>());
      if (!cache.unl[ev->ml].count(ev->thread_id()))
        cache.unl[ev->ml].emplace(ev->thread_id(), std::vector<const ZEvent *>());
      cache.unl[ev->ml][ev->thread_id()].push_back(ev);
    }

  } // end of loop for traversing trace and creating nodes

  shrink();
  po.shrink();

  #ifndef NDEBUG
  assert(size() == cur_evidx.size());
  for (auto& thraux_numevents : cur_evidx) {
    assert(thraux_numevents.second == (*this)(thraux_numevents.first).size()
           && "Didn't go through entire original part of the (*this)");
  }
  #endif

  // EDGES - spawns
  for (const ZEvent *spwn : spawns) {
    if (!proc_seq_to_thread_id.count(spwn->childs_cpid.get_proc_seq()))
      continue;
    unsigned child_thr = proc_seq_to_thread_id.at(spwn->childs_cpid.get_proc_seq());
    for (int aux : auxes(child_thr)) {
      assert(hasThreadAux(child_thr, aux));
      assert(!(*this)(child_thr, aux).empty());
      const ZEvent *nthr = getEvent(child_thr, aux, 0);

      assert(!po.hasEdge(nthr, spwn));
      if (!po.hasEdge(spwn, nthr))
        po.addEdge(spwn, nthr);
    }
  }

  // EDGES - joins
  for (const ZEvent *jn : joins) {
    if (!proc_seq_to_thread_id.count(jn->childs_cpid.get_proc_seq()))
      continue;
    unsigned child_thr = proc_seq_to_thread_id.at(jn->childs_cpid.get_proc_seq());
    assert(child_thr < proc_seq_to_thread_id.size());
    for (int aux : auxes(child_thr)) {
      assert(hasThreadAux(child_thr, aux));
      assert(!(*this)(child_thr, aux).empty());
      unsigned lastev_idx = (*this)(child_thr, aux).size() - 1;
      const ZEvent *wthr = getEvent(child_thr, aux, lastev_idx);

      assert(!po.hasEdge(jn, wthr));
      if (!po.hasEdge(wthr, jn))
        po.addEdge(wthr, jn);
    }
  }

  // EDGES - mutex inits
  for (auto& in : mutex_inits) {
    const SymAddrSize& loc_init = in.first;
    const ZEvent *ev_init = in.second;

    for (auto& tid_ev_first : mutex_first[loc_init]) {
      const ZEvent *ev_first = tid_ev_first.second;

      assert(!po.hasEdge(ev_first, ev_init));
      if (!po.hasEdge(ev_init, ev_first))
        po.addEdge(ev_init, ev_first);
    }
  }

  // EDGES - mutex destroys
  for (auto& de : mutex_destroys) {
    const SymAddrSize& loc_destroy = de.first;
    const ZEvent *ev_destroy = de.second;

    for (auto& tid_ev_last : mutex_last[loc_destroy]) {
      const ZEvent *ev_last = tid_ev_last.second;

      assert(!po.hasEdge(ev_destroy, ev_last));
      if (!po.hasEdge(ev_last, ev_destroy))
        po.addEdge(ev_last, ev_destroy);
    }
  }
  end_err("ending trace-to-po");
}


void ZGraph::add_reads_from_edges(const ZAnnotation& annotation)
{
  start_err("starting reads-from-edges");
  // Reads
  for (auto it = annotation.read_begin(); it != annotation.read_end(); ++it) {
    assert(hasEvent(it->first));
    const ZEvent * const read = getEvent(it->first);
    assert(isRead(read) && !it->second.cpid().is_auxiliary());
    if (it->second != ZEventID(true) && read->cpid() != it->second.cpid()) {
      assert(hasEvent(it->second));
      const ZEvent * const wM = getEvent(it->second)->write_other_ptr;
      assert(isWriteM(wM) && sameMl(read, wM));
      assert(!po.hasEdge(read, wM) && "Inconsistent annotation");
      if (!po.hasEdge(wM, read)) po.addEdge(wM, read);
    }
  }
  // Locks
  for (auto it = annotation.lock_begin(); it != annotation.lock_end(); ++it) {
    assert(hasEvent(it->first));
    const ZEvent * const lock = getEvent(it->first);
    assert(isLock(lock) && !it->second.cpid().is_auxiliary());
    if (it->second != ZEventID(true) && lock->cpid() != it->second.cpid()) {
      assert(hasEvent(it->second));
      const ZEvent * const unlock = getEvent(it->second);
      assert(isUnlock(unlock) && sameMl(lock, unlock));
      assert(!po.hasEdge(lock, unlock) && "Inconsistent annotation");
      if (!po.hasEdge(unlock, lock)) po.addEdge(unlock, lock);
    }
  }
  end_err("ending reads-from-edges");
}


/* *************************** */
/* MAIN ALGORITHM              */
/* *************************** */

std::set<ZEventID> ZGraph::get_mutations
(const ZEvent * const ev, const ZEventID base_obs,
 int mutations_only_from_idx) const
{
  start_err("starting get-mutations");
  assert(isRead(ev) || isLock(ev));
  assert(hasEvent(ev));
  std::set<ZEventID> res;
  const ZEvent * const ev_before = (ev->event_id() > 0)
    ? getEvent(ev->thread_id(), ev->aux_id(), ev->event_id() - 1) : nullptr;
  const ZEvent * spawn = proc_seq_to_spawn.count(ev->cpid().get_proc_seq()) ?
  proc_seq_to_spawn.at(ev->cpid().get_proc_seq()) : nullptr;
  assert(spawn || ev->cpid().get_proc_seq() == CPid().get_proc_seq());
  assert(!spawn || po.hasEdge(spawn, ev));

  std::unordered_set<const ZEvent *> notCovered;
  std::unordered_set<const ZEvent *> mayBeCovered;
  // From the value (evid) and below, everything is covered
  // from read/lock by some other write/unlock
  std::unordered_map<std::pair<unsigned, int>, int> covered;
  for (unsigned covthr : get_threads())
    for (int covaux : auxes(covthr))
      covered.emplace(std::pair<unsigned, int>(covthr, covaux), -1);

  // Handle other threads
  for (unsigned tid : get_threads()) {
    int auxid = isLock(ev) ? -1 : auxForMl(ev->ml, tid);
    if (ev->thread_id() == tid || (auxid < 0 && isRead(ev)))
      continue;

    if (isLock(ev)) {
      if (!getCache().unl.count(ev->ml) ||
          !getCache().unl.at(ev->ml).count(tid))
        continue;
      const ZEvent * last_unlock_before_lock = nullptr;
      for (const ZEvent * unlock : getCache().unl.at(ev->ml).at(tid)) {
        assert(isUnlock(unlock) && tid == unlock->thread_id() && sameMl(ev, unlock));
        if (po.hasEdge(ev, unlock))
          break;
        else if (ev_before && po.hasEdge(unlock, ev_before))
          last_unlock_before_lock = unlock;
        else {
          assert(!notCovered.count(unlock));
          notCovered.emplace(unlock);
        }
      }
      if (last_unlock_before_lock) {
        assert(!mayBeCovered.count(last_unlock_before_lock));
        mayBeCovered.emplace(last_unlock_before_lock);
        // Update cover
        for (unsigned covthr : get_threads()) {
          if (covthr != last_unlock_before_lock->thread_id()) {
            for (int covaux : auxes(covthr)) {
              int newcov = po.pred(last_unlock_before_lock, covthr, covaux).second;
              auto covta = std::pair<unsigned, int>(covthr, covaux);
              assert(covered.count(covta));
              if (newcov > covered[covta])
                covered[covta] = newcov;
            }
          }
        }
      }
    } else {
      assert(isRead(ev) && auxid >= 0);
      int su = po.succ(ev, tid, auxid).second;
      if (su >= (int) (*this)(tid, auxid).size()) {
        assert(su == INT_MAX);
        su = (int) (*this)(tid, auxid).size();
      }
      // Dumb
      /*
      int curr = su - 1;
      while (curr >= 0) {
        const ZEvent *rem = (*this)(tid, auxid)[curr];
        assert(isWriteM(rem) && !po.hasEdge(ev, rem));
        if (po.hasEdge(rem, ev)) {
          if (sameMl(ev, rem)) {
            mayBeCovered.emplace(rem);
            // Update cover
            for (unsigned covthr : get_threads()) {
              if (covthr != rem->thread_id()) {
                for (int covaux : auxes(covthr)) {
                  int newcov = po.pred(rem, covthr, covaux).second;
                  auto covta = std::pair<unsigned, int>(covthr, covaux);
                  assert(covered.count(covta));
                  if (newcov > covered[covta])
                    covered[covta] = newcov;
                }
              }
            }
            // Break
            break;
          }
          curr--;
        } else {
          if (sameMl(ev, rem))
            notCovered.emplace(rem);
          curr--;
        }
      }
      */
      // TailW index to cache.wm
      int wm_idx = getTailWindex(ev->ml, tid, su - 1);
      assert(wm_idx == getLatestNotAfterIndex(ev, tid, po));
      while (wm_idx > -1) {
        const ZEvent *rem = cache.wm.at(ev->ml).at(tid).at(wm_idx);
        assert(rem && isWriteM(rem) && sameMl(ev, rem));
        assert(!po.hasEdge(ev, rem));
        if ((ev_before && po.hasEdge(rem, ev_before)) ||
            (spawn && po.hasEdge(rem, spawn))) {
          // Others in this thread are covered from ev by rem
          assert((spawn && po.hasEdge(rem, spawn)) ||
                 (ev_before && (int) rem->event_id() <= po.pred(ev_before, tid, auxid).second));
          mayBeCovered.emplace(rem);
          // Update cover
          for (unsigned covthr : get_threads()) {
            if (covthr != rem->thread_id()) {
              for (int covaux : auxes(covthr)) {
                int newcov = po.pred(rem, covthr, covaux).second;
                auto covta = std::pair<unsigned, int>(covthr, covaux);
                assert(covered.count(covta));
                if (newcov > covered[covta])
                  covered[covta] = newcov;
              }
            }
          }
          // Break
          break;
        } else {
          notCovered.emplace(rem);
        }
        --wm_idx;
      }
    }
  }

  // Handle thread of read/lock
  if (isLock(ev)) {
    const ZEvent * last_unlock_before_lock = nullptr;
    if (getCache().unl.count(ev->ml) && getCache().unl.at(ev->ml).count(ev->thread_id())) {
      for (const ZEvent * unlock : getCache().unl.at(ev->ml).at(ev->thread_id())) {
        assert(isUnlock(unlock) && ev->thread_id() == unlock->thread_id() &&
               sameMl(ev, unlock) && !unlock->cpid().is_auxiliary());
        if (unlock->event_id() > ev->event_id())
          break;
        else {
          assert(unlock->event_id() < ev->event_id());
          assert(!last_unlock_before_lock ||
                 last_unlock_before_lock->event_id() < unlock->event_id());
          last_unlock_before_lock = unlock;
        }
      }
    }
    if (last_unlock_before_lock) {
      // Update cover
      for (unsigned covthr : get_threads()) {
        if (covthr != last_unlock_before_lock->thread_id()) {
          for (int covaux : auxes(covthr)) {
            int newcov = po.pred(last_unlock_before_lock, covthr, covaux).second;
            auto covta = std::pair<unsigned, int>(covthr, covaux);
            assert(covered.count(covta));
            if (newcov > covered[covta])
              covered[covta] = newcov;
          }
        }
      }

      // Check if not covered by some remote
      bool localCovered = false;
      for (const auto& rem : mayBeCovered) {
        assert(isUnlock(rem) && ev_before && po.hasEdge(rem, ev_before));
        if (po.hasEdge(last_unlock_before_lock, rem)) {
          localCovered = true;
          break;
        }
      }
      // Add obs if not covered, not base_obs, from_idx
      if (!localCovered && base_obs != last_unlock_before_lock->id() &&
          last_unlock_before_lock->trace_id() >= mutations_only_from_idx) {
        assert(!res.count(last_unlock_before_lock->id()));
        res.emplace(last_unlock_before_lock->id());
      }
    } else {
       // No local unlock for lock
      if (mayBeCovered.empty() && base_obs != ZEventID(true) &&
          mutations_only_from_idx == -1) {
        // Consider initial event
        assert(!res.count(ZEventID(true)));
        res.emplace(ZEventID(true)); // initial
      }
    }
  } else {
    assert(isRead(ev));
    const ZEvent *localB = cache.readWB.at(ev);
    if (localB) {
      // There is a local write for read
      assert(isWriteB(localB) && sameMl(localB, ev) &&
             localB->thread_id() == ev->thread_id() &&
             localB->aux_id() == -1 &&
             localB->event_id() < ev->event_id());
      const ZEvent *localM = localB->write_other_ptr;
      // Update cover caused by localB
      for (unsigned covthr : get_threads()) {
        if (covthr != localB->thread_id()) {
          for (int covaux : auxes(covthr)) {
            int newcov = po.pred(localB, covthr, covaux).second;
            auto covta = std::pair<unsigned, int>(covthr, covaux);
            assert(covered.count(covta));
            if (newcov > covered[covta])
              covered[covta] = newcov;
          }
        }
      }
      // Delete remotes that happen *before* localM
      // (no chance for them, those after localM are safe)
      assert(isWriteM(localM) && sameMl(localB, localM) &&
             po.hasEdge(localB, localM));
      for (auto it = notCovered.begin(); it != notCovered.end(); ) {
        const ZEvent *rem = *it;
        assert(isWriteM(rem));
        assert(!ev_before || !po.hasEdge(rem, ev_before));
        if (po.hasEdge(rem, localM))
          it = notCovered.erase(it);
        else
          ++it;
      }
      for (auto it = mayBeCovered.begin(); it != mayBeCovered.end(); ) {
        const ZEvent *rem = *it;
        assert(isWriteM(rem));
        assert((ev_before && po.hasEdge(rem, ev_before)) ||
               (spawn && po.hasEdge(rem, spawn)));
        if (po.hasEdge(rem, localM))
          it = mayBeCovered.erase(it);
        else
          ++it;
      }
      // Check if not covered by some remote
      bool localCovered = false;
      for (const auto& rem : mayBeCovered) {
        assert(isWriteM(rem));
        assert((ev_before && po.hasEdge(rem, ev_before)) ||
               (spawn && po.hasEdge(rem, spawn)));
        if (po.hasEdge(localM, rem)) {
          localCovered = true;
          break;
        }
      }
      // Add obs if not covered, not base_obs, from_idx
      if (!localCovered && base_obs != localB->id() &&
          localB->trace_id() >= mutations_only_from_idx) {
        assert(!res.count(localB->id()));
        res.emplace(localB->id());
      }
    } else {
      // No local write for read
      if (mayBeCovered.empty() && base_obs != ZEventID(true) &&
          mutations_only_from_idx == -1) {
        // Consider initial event
        assert(!res.count(ZEventID(true)));
        res.emplace(ZEventID(true)); // initial
      }
    }
  }


  // Take candidates unordered with read/lock
  for (const auto& not_cov : notCovered) {
    const ZEvent * obs_ev = not_cov;
    assert(isWriteM(obs_ev) || isUnlock(obs_ev));
    if (isWriteM(obs_ev)) {
      obs_ev = obs_ev->write_other_ptr;
      assert(isWriteB(obs_ev));
    }
    assert(sameMl(obs_ev, ev));
    // Add if not base_obs, from_idx
    if (base_obs != obs_ev->id() &&
        obs_ev->trace_id() >= mutations_only_from_idx) {
      assert(!res.count(obs_ev->id()));
      res.emplace(obs_ev->id());
    }
  }

  // Take candidates that happen before read/lock
  for (const auto& may_be : mayBeCovered) {
    const ZEvent * obs_ev = may_be;
    assert(isWriteM(obs_ev) || isUnlock(obs_ev));
    auto covta = std::pair<unsigned, int>
      (obs_ev->thread_id(), obs_ev->aux_id());
    assert(covered.count(covta));
    if ((int) obs_ev->event_id() > covered[covta]) {
      // Not covered
      if (isWriteM(obs_ev)) {
        obs_ev = obs_ev->write_other_ptr;
        assert(isWriteB(obs_ev));
      }
      assert(sameMl(obs_ev, ev));
      // Add if not base_obs, from_idx
      if (base_obs != obs_ev->id() &&
          obs_ev->trace_id() >= mutations_only_from_idx) {
        assert(!res.count(obs_ev->id()));
        res.emplace(obs_ev->id());
      }
    }
  }
  end_err("ending get-mutations");
  return res;
}


std::pair<std::set<int>, std::set<const ZEvent *>>
ZGraph::get_causes_after
(const ZEvent * const readlock, const ZEventID& mutation,
 const std::vector<std::unique_ptr<ZEvent>>& tau) const
{
  assert(isRead(readlock) || isLock(readlock));
  if (mutation == ZEventID(true))
    return {std::set<int>(), std::set<const ZEvent *>()};

  std::set<int> causes_all_idx;
  std::set<const ZEvent *> causes_readslocks;

  assert(hasEvent(mutation));

  const ZEvent * obs = getEvent(mutation);
  assert(isWriteB(obs) || isUnlock(obs));
  if (isWriteB(obs)) {
    obs = obs->write_other_ptr;
    assert(isWriteM(obs));
  }
  assert(sameMl(obs, readlock));

  for (int i = readlock->trace_id() + 1; i <= obs->trace_id(); ++i) {
    assert(i >= 1 && i < tau.size());
    const ZEvent * const ev = tau[i].get();
    assert(hasEvent(ev));
    if (*ev == *obs || getPo().hasEdge(ev, obs)) {
      assert(!causes_all_idx.count(ev->trace_id()));
      causes_all_idx.emplace(ev->trace_id());
      if (isRead(ev) || isLock(ev)) {
        assert(!causes_readslocks.count(ev));
        causes_readslocks.emplace(ev);
      }
    }
  }
  return {causes_all_idx, causes_readslocks};
}


bool ZGraph::is_causal_readlock_or_mutation
(const ZEvent * const causal, const ZEvent * const readlock,
 const ZEvent * const mut_ev) const
{
  assert(isRead(causal) || isLock(causal));
  assert(isRead(readlock) || isLock(readlock));
  assert(hasEvent(causal) && hasEvent(readlock));
  assert(!mut_ev || isWriteB(mut_ev) || isUnlock(mut_ev));
  assert(!mut_ev || hasEvent(mut_ev));
  if (*causal == *readlock) {
    // We return false for readlock itself
    // (it is added separately in explorer)
    return false;
  }
  if (mut_ev && po.hasEdge(causal, mut_ev))
    return true; // Causal past of mutation
  // Check causal past of readlock
  const ZEvent * const ev_before = (readlock->event_id() > 0)
    ? getEvent(readlock->thread_id(), readlock->aux_id(), readlock->event_id() - 1) : nullptr;
  const ZEvent * spawn = proc_seq_to_spawn.count(readlock->cpid().get_proc_seq()) ?
  proc_seq_to_spawn.at(readlock->cpid().get_proc_seq()) : nullptr;
  assert(spawn || readlock->cpid().get_proc_seq() == CPid().get_proc_seq());
  assert(!spawn || po.hasEdge(spawn, readlock));
  if (ev_before && (*causal == *ev_before || po.hasEdge(causal, ev_before)))
    return true;
  if (spawn && po.hasEdge(causal, spawn))
    return true;
  return false;
}


int ZGraph::getTailWindex(const SymAddrSize& ml, unsigned thr, int evX) const
{
  assert(threads_auxes.count(thr));
  assert(cache.wm.count(ml));
  if (evX < 0)
    return -1;
  unsigned ev = evX;
  if (!cache.wm.at(ml).count(thr))
    return -1;

  const auto& writes = cache.wm.at(ml).at(thr);
  if (writes.empty())
    return -1;

  int low = 0;
  if (ev < writes[low]->event_id())
    return -1;

  int high = writes.size() - 1;
  if (ev >= writes[high]->event_id())
    return high;

  assert(low < high);
  if (low + 1 == high) {
    assert(ev >= writes[low]->event_id() &&
           ev < writes[high]->event_id());
    return low;
  }

  // Low represents something that
  // can possibly be the answer
  // High represents something that
  // is above the answer
  // Do binary search
  while (true) {
    assert(low + 1 < high);
    assert(ev >= writes[low]->event_id() &&
           ev < writes[high]->event_id());
    int mid = ((high - low) / 2) + low;
    assert(low < mid && mid < high);

    if (ev >= writes[mid]->event_id())
      low = mid;
    else
      high = mid;

    if (low + 1 == high) {
      assert(ev >= writes[low]->event_id() &&
             ev < writes[high]->event_id());
      return low;
    }
  }
}


const ZEvent * const ZGraph::getTailW
(const SymAddrSize& ml, unsigned thr, int evX) const
{
  int idx = getTailWindex(ml, thr, evX);
  if (idx == -1)
    return nullptr;

  assert(idx < (int) cache.wm.at(ml).at(thr).size());
  auto res = cache.wm.at(ml).at(thr)[idx];
  assert(isWriteM(res) && res->ml == ml &&
         res->thread_id() == thr && res->aux_id() != -1 &&
         evX >= 0 && res->event_id() <= (unsigned) evX);
  return res;
}


int ZGraph::getLatestNotAfterIndex(const ZEvent *read, unsigned thr, const ZPartialOrder& partial) const
{
  assert(isRead(read));
  assert(hasEvent(read));
  assert(threads_auxes.count(thr));
  assert(thr != read->thread_id());
  assert(cache.wm.count(read->ml));

  int aux = auxForMl(read->ml, thr);
  if (aux == -1) {
    assert(!cache.wm.at(read->ml).count(thr));
    return -1;
  }
  assert(aux == 0 || cache.wm.at(read->ml).count(thr));

  int su = partial.succ(read, thr, aux).second;
  if (su >= (int) (*this)(thr, aux).size()) {
    assert(su == INT_MAX);
    su = (int) (*this)(thr, aux).size();
  }
  // TailW index to cache.wm
  return getTailWindex(read->ml, thr, su - 1);
}


const ZEvent * const ZGraph::getLatestNotAfter(const ZEvent *read, unsigned thr, const ZPartialOrder& partial) const
{
  assert(isRead(read));
  assert(hasEvent(read));
  int idx = getLatestNotAfterIndex(read, thr, partial);
  if (idx == -1)
    return nullptr;

  assert(idx < (int) cache.wm.at(read->ml).at(thr).size());
  auto res = cache.wm.at(read->ml).at(thr)[idx];
  assert(isWriteM(res) && sameMl(res, read) &&
         res->thread_id() == thr && res->aux_id() != -1 &&
         !partial.hasEdge(read, res));
  return res;
}


const ZEvent * const ZGraph::getLocalBufferW(const ZEvent *read) const
{
  assert(isRead(read));
  assert(hasEvent(read));
  assert(cache.readWB.count(read));
  const ZEvent *local = cache.readWB.at(read);
  assert(!local || isWriteB(local));
  return local;
}
