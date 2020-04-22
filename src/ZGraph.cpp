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
#include "ZBuilderSC.h"

static const bool DEBUG = false;
#include "ZDebug.h"


// Empty
ZGraph::ZGraph()
  : tso(true),
    init(ZEvent(true)),
    lines(),
    thread_aux_to_line_id(),
    threads_auxes(),
    pso_thr_mlauxes(),
    proc_seq_to_thread_id(),
    event_to_position(),
    po(*this),
    cache()
{
  assert(empty());
}


// Initial
ZGraph::ZGraph(const std::vector<ZEvent>& trace, bool tso)
  : tso(tso),
    init(ZEvent(true)),
    lines(),
    thread_aux_to_line_id(),
    threads_auxes(),
    pso_thr_mlauxes(),
    proc_seq_to_thread_id(),
    event_to_position(),
    po(*this),
    cache()
{
  assert(empty());
  traceToPO(trace, nullptr);
}


// Moving
ZGraph::ZGraph(ZGraph&& oth)
  : tso(oth.tso),
    init(std::move(oth.init)),
    lines(std::move(oth.lines)),
    thread_aux_to_line_id(std::move(oth.thread_aux_to_line_id)),
    threads_auxes(std::move(oth.threads_auxes)),
    pso_thr_mlauxes(std::move(oth.pso_thr_mlauxes)),
    proc_seq_to_thread_id(std::move(oth.proc_seq_to_thread_id)),
    event_to_position(std::move(oth.event_to_position)),
    po(std::move(oth.po)),
    cache(std::move(oth.cache))
{
  assert(oth.empty());
}


// Extending
// Partial order that will be moved as original
// Trace and annotation that will extend this copy of the graph
ZGraph::ZGraph
(const ZGraph& oth,
 ZPartialOrder&& po,
 const std::vector<ZEvent>& trace,
 const ZAnnotation& annotation)
  : tso(oth.tso),
    init(ZEvent(true)),
    lines(oth.lines),
    thread_aux_to_line_id(oth.thread_aux_to_line_id),
    threads_auxes(oth.threads_auxes),
    pso_thr_mlauxes(oth.pso_thr_mlauxes),
    proc_seq_to_thread_id(oth.proc_seq_to_thread_id),
    event_to_position(oth.event_to_position),
    po(std::move(po), *this),
    cache()
{
  traceToPO(trace, &annotation);
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


const ZEvent * ZGraph::event(unsigned thread_id, int aux_id, unsigned event_id) const
{
  assert(thread_id != INT_MAX && "Called for initial event");
  assert(hasThreadAux(thread_id, aux_id));
  const LineT& line = this->operator()(thread_id, aux_id);
  assert(event_id < line.size());
  assert(has_event(line[event_id]));
  return line[event_id];
}


const ZEvent * ZGraph::event(const ZObs& obs) const
{
  return event(obs.thr, -1, obs.ev);
}


const ZEvent * ZGraph::last_unlock_of_mutex(const ZObs& obs) const
{
  unsigned curEv = obs.ev;
  auto lock = event(obs.thr, -1, curEv);
  assert(is_lock(lock));

  while (curEv + 1 < (*this)(obs.thr, -1).size()) {
    curEv++;
    auto res = event(obs.thr, -1, curEv);
    if (is_unlock(res) && same_ml(lock, res))
      return res;
  }

  return nullptr;
}


void ZGraph::add_line(const ZEvent *ev)
{
  assert(ev);
  assert(!has_event(ev));
  assert(ev->thread_id() < 20 && "Thread ID not set up yet");
  assert(!hasThreadAux(ev->thread_id(), ev->aux_id()));

  auto key = std::pair<unsigned, int>(ev->thread_id(), ev->aux_id());
  thread_aux_to_line_id.emplace(key, lines.size());
  assert(ev->thread_id() <= threads_auxes.size());
  if (ev->thread_id() == threads_auxes.size())
    threads_auxes.push_back(std::set<int>());
  assert(!threads_auxes.at(ev->thread_id()).count(ev->aux_id()));
  threads_auxes[ev->thread_id()].emplace(ev->aux_id());
  lines.push_back(LineT());
  lines.back().reserve(16);
}


void ZGraph::add_event(const ZEvent *ev)
{
  assert(ev);
  assert(!has_event(ev));
  assert(ev->thread_id() < 20 && "Thread ID not set up yet");
  assert(hasThreadAux(ev->thread_id(), ev->aux_id()));

  auto key = std::pair<unsigned, int>(ev->thread_id(), ev->aux_id());
  auto it = thread_aux_to_line_id.find(key);
  assert(it != thread_aux_to_line_id.end());
  LineT& line = lines[it->second];

  unsigned event_id = line.size();
  assert(event_id == ev->event_id());
  assert(tso || ev->aux_id() == -1 || line.empty() || same_ml(line[0], ev));
  event_to_position.emplace(ev, std::pair<unsigned,unsigned>(it->second, event_id));
  line.push_back(ev);
  assert(has_event(ev));
}


void ZGraph::replace_event(const ZEvent *oldEv, const ZEvent *newEv)
{
  assert(oldEv && newEv);
  assert(oldEv != newEv && (*oldEv == *newEv)
         && "different pointers but same thread/aux/eventID");
  assert(oldEv->kind == newEv->kind &&
         oldEv->cpid() == newEv->cpid() &&
         oldEv->ml == newEv->ml);
  assert(has_event(oldEv) && !has_event(newEv));

  auto line_event = event_to_position[oldEv];
  event_to_position.erase(oldEv);
  event_to_position.emplace(newEv, line_event);

  assert(line_event.first < lines.size() &&
         line_event.second < lines[line_event.first].size() &&
         lines[line_event.first][line_event.second] == oldEv);
  lines[line_event.first][line_event.second] = newEv;
  assert(!has_event(oldEv) && has_event(newEv));
}


void ZGraph::shrink()
{
  lines.shrink_to_fit();
  for (unsigned i=0; i<lines.size(); ++i)
    lines[i].shrink_to_fit();
}


unsigned ZGraph::line_id(unsigned thread_id, int aux_id) const
{
  assert(thread_id != INT_MAX && "Called for initial event");
  assert(hasThreadAux(thread_id, aux_id));
  auto key = std::pair<unsigned,int>(thread_id, aux_id);
  auto it = thread_aux_to_line_id.find(key);
  assert(it != thread_aux_to_line_id.end());
  return it->second;
}


unsigned ZGraph::line_id(const ZEvent *ev) const
{
  assert(ev);
  assert(!is_initial(ev) && "Called for initial event");
  return line_id(ev->thread_id(), ev->aux_id());
}


bool ZGraph::hasThreadAux(std::pair<unsigned, int> ids) const
{
  assert(ids.first != INT_MAX && "Called for initial event");
  assert(threads_auxes.size() <= ids.first ||
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


std::vector<unsigned> ZGraph::real_sizes_minus_one() const
{
  std::vector<unsigned> res;
  for (unsigned thr = 0; thr < number_of_threads(); ++thr) {
    assert(!(*this)(thr, -1).empty());
    res.push_back((*this)(thr, -1).size() - 1);
  }
  return res;
}


unsigned ZGraph::number_of_threads() const
{
  return threads_auxes.size();
}


const std::set<int>& ZGraph::auxes(unsigned thread_id) const
{
  assert(thread_id != INT_MAX && "Called for initial event");
  assert(thread_id < threads_auxes.size());
  return threads_auxes[thread_id];
}


int ZGraph::auxForMl(const SymAddrSize& ml, unsigned thr) const
{
  assert(thr < number_of_threads());
  auto axs = auxes(thr);
  assert(!axs.empty());
  if (axs.size() == 1) {
    // No write events in this thread
    return -1;
  }
  if (axs.size() == 2 && tso) {
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
  assert(!tso);
  assert(pso_thr_mlauxes.count(thr));
  for (auto aux : axs) {
    if (aux != -1) {
      assert(!((*this)(thr, aux).empty()));
      const ZEvent *first = (*this)(thr, aux)[0];
      assert(is_writeM(first));
      if (first->ml == ml)
        return aux;
    }
  }
  // This thread has no event for this ml
  return -1;
}


int ZGraph::psoGetAux(const ZEvent* writeM)
{
  assert(!tso && is_writeM(writeM));
  unsigned thr_idx = writeM->thread_id();
  assert(thr_idx < 100);
  int aux_idx = -1;
  if (!pso_thr_mlauxes.count(thr_idx)) {
    // First aux
    pso_thr_mlauxes.emplace
      (thr_idx, std::vector<SymAddrSize>());
    pso_thr_mlauxes[thr_idx].push_back(writeM->ml);
    aux_idx = 0;
  } else {
    for (unsigned ax = 0; ax < pso_thr_mlauxes[thr_idx].size(); ++ax)
      if (pso_thr_mlauxes[thr_idx][ax] == writeM->ml) {
        aux_idx = ax;
        break;
      }
    if (aux_idx == -1) {
      // New aux
      aux_idx = pso_thr_mlauxes[thr_idx].size();
      pso_thr_mlauxes[thr_idx].push_back(writeM->ml);
    }
  }
  assert(aux_idx >= 0);
  return aux_idx;
}


// <threadID, added_with_this_call?>
std::pair<unsigned, bool> ZGraph::getThreadID(const std::vector<int>& proc_seq)
{
  auto it = proc_seq_to_thread_id.find(proc_seq);
  if (it == proc_seq_to_thread_id.end()) {
    unsigned res = proc_seq_to_thread_id.size();
    proc_seq_to_thread_id.emplace_hint(it, proc_seq, res);
    return {res, true};
  };
  return {it->second, false};
}


// <threadID, added_with_this_call?>
std::pair<unsigned, bool> ZGraph::getThreadID(const ZEvent *ev)
{
  assert(ev);
  assert(!is_initial(ev) && "Called for initial event");
  return getThreadID(ev->cpid().get_proc_seq());
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
  assert(!is_initial(ev) && "Called for initial event");
  return getThreadIDnoAdd(ev->cpid().get_proc_seq());
}


bool ZGraph::has_event(const ZEvent *ev) const
{
  assert(ev);
  if (ev == initial())
    return true;
  auto it = event_to_position.find(ev);
  if (it == event_to_position.end())
    return false;
  assert(lines[it->second.first][it->second.second] == ev);
  assert(line_id(ev) == it->second.first);
  assert(ev->event_id() == it->second.second);
  return true;
}


/* *************************** */
/* GRAPH EXTENSION             */
/* *************************** */

// Extends this graph so it corresponds to 'trace'
// Check the header file for the method description
void ZGraph::traceToPO(const std::vector<ZEvent>& trace,
                       const ZAnnotation *annotationPtr)
{
  err_msg("starting trace-to-po");
  assert(!trace.empty());
  assert(trace.size() >= events_size());

  std::unordered_map<int, std::pair<unsigned, int>> ipid_to_thraux;
  ipid_to_thraux.reserve(8);

  std::unordered_map<std::pair<unsigned, int>, unsigned>
    cur_evidx;
  cur_evidx.reserve(8);

  std::unordered_set<unsigned> has_unannotated_read_or_lock;
  has_unannotated_read_or_lock.reserve(8);

  std::unordered_set<SymAddrSize> found_last_lock_for_location;
  found_last_lock_for_location.reserve(8);

  std::unordered_set<unsigned> forbidden_processes_ipid;
  forbidden_processes_ipid.reserve(8);

  std::unordered_set<std::vector<int>> proc_seq_within_po;
  proc_seq_within_po.reserve(8);
  proc_seq_within_po.insert(trace[0].cpid().get_proc_seq());

  std::unordered_set<std::pair<unsigned, int>> thraux_with_event_we_dont_add;
  thraux_with_event_we_dont_add.reserve(8);

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

  for (auto traceit = trace.begin(); traceit != trace.end(); ++traceit) {
    const ZEvent *ev = &(*traceit);

    if (forbidden_processes_ipid.count(ev->iid.get_pid())) {
        // This is a forbidden process because it
        // happens after a so-far unannotated node
        continue;
    }
    if (!proc_seq_within_po.count(ev->cpid().get_proc_seq())) {
      // This is a forbidden process because it
      // is created after a so-far unannotated node
      forbidden_processes_ipid.insert(ev->iid.get_pid());
      continue;
    }

    unsigned thr_idx = INT_MAX;

    // Check if this process is already known
    auto ipidit = ipid_to_thraux.find(ev->iid.get_pid());
    if (tso && ipidit != ipid_to_thraux.end()) {
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
    if (!tso && is_writeM(ev)) {
      // Handle aux for pso
      int aux_idx = psoGetAux(ev);
      ev->_aux_id = aux_idx;
    }

    auto thraux = std::pair<unsigned, int>(ev->thread_id(), ev->aux_id());

    // Check possible cases why we cannot add this event
    //
    // Check if we already haven't added a thraux-predecessor of this event
    if (thraux_with_event_we_dont_add.count(thraux)) {
      assert(!hasThreadAux(thraux) ||
             cur_evidx.at(thraux) == (*this)(thraux).size());
      continue;
    }
    // Check if thraux already has an unannotated read
    if (has_unannotated_read_or_lock.count(thraux.first) &&
        thraux.second == -1) {
      // This event happens after something we need to
      // annotate first, so we don't add it into the graph
      assert(!hasThreadAux(thraux) ||
             cur_evidx.at(thraux) == (*this)(thraux).size());
      thraux_with_event_we_dont_add.insert(thraux);
      continue;
    }
    // Joins of threads with an event we didn't include also
    // can not be included and neither any of their successors
    if (is_join(ev)) {
      bool dontInclude = false;
      auto childs_proc_it = proc_seq_within_po.find(ev->childs_cpid.get_proc_seq());
      if (childs_proc_it == proc_seq_within_po.end()) {
        dontInclude = true;
      } else {
        auto thr_added = getThreadID(ev->childs_cpid.get_proc_seq());
        assert(!thr_added.second);
        auto childs_thraux = std::pair<unsigned, int>(thr_added.first, -1);
        // It is enough to check the presence of the real thread
        if (thraux_with_event_we_dont_add.count(childs_thraux))
          dontInclude = true;
      }
      if (dontInclude) {
        assert(!hasThreadAux(thraux) ||
               cur_evidx.at(thraux) == (*this)(thraux).size());
        thraux_with_event_we_dont_add.insert(thraux);
        continue;
      }
    }
    // Memory-writes cannot be added if their buffer-write
    // was not added
    if (is_writeM(ev)) {
      assert(ev->aux_id() != -1 && "Only auxiliary threads");
      if (!store_buffer.count(ev->thread_id()) ||
          !store_buffer[ev->thread_id()].count(ev->ml) ||
          store_buffer[ev->thread_id()][ev->ml].empty()) {
        assert(!hasThreadAux(thraux) ||
               cur_evidx.at(thraux) == (*this)(thraux).size());
        thraux_with_event_we_dont_add.insert(thraux);
        continue;
      }
    }

    // We will add the event
    // First, add line if not already present
    if (!hasThreadAux(thraux)) {
      // A new thread/aux appeared,
      // we've checked above that it is allowed
      assert(!cur_evidx.count(thraux));
      cur_evidx.emplace(thraux, 0);

      add_line(ev);
      po.add_line(ev);
    }

    if (!cur_evidx.count(thraux)) {
      assert(ev->event_id() == 0);
      cur_evidx.emplace(thraux, 0);
    }
    unsigned ev_idx = cur_evidx[thraux];

    assert(ev->event_id() == ev_idx);
    if (ev_idx == (*this)(thraux).size()) {
      // New event
      add_event(ev);
      po.add_event(ev);

      // XXX: function pointer calls not handled
      if (is_spawn(ev))
        spawns.push_back(ev);
      if (is_join(ev))
        joins.push_back(ev);

    } else {
      // Already known event
      const ZEvent *oldev = event(ev->thread_id(),
                                     ev->aux_id(),
                                     ev->event_id());
      replace_event(oldev, ev);
    }

    ++cur_evidx[thraux];

    // Handle Fence
    //
    if (ev->fence) {
      if (store_buffer.count(ev->thread_id()))
        for (const auto& ml_last : last_mwrite[ev->thread_id()]) {
          auto last = ml_last.second;
          assert(is_writeM(last));
          assert(!po.has_edge(ev, last));
          if (!po.has_edge(last, ev))
            po.add_edge(last, ev);
        }
    }

    // Handle Spawn
    //
    if (is_spawn(ev)) {
      proc_seq_within_po.insert(ev->childs_cpid.get_proc_seq());
      assert(number_of_threads() > 0);
      assert(ev->fence);
    }

    // Handle Buffer-Write
    //
    if (is_writeB(ev)) {
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
    //
    if (is_writeM(ev)) {
      assert(ev->aux_id() != -1 && "Only auxiliary threads");
      // Store buffer
      assert(store_buffer.count(ev->thread_id()) &&
             store_buffer[ev->thread_id()].count(ev->ml));
      assert(!ev->write_other_ptr);
      ev->write_other_ptr = store_buffer[ev->thread_id()][ev->ml].front();
      assert(is_writeB(ev->write_other_ptr) && same_ml(ev, ev->write_other_ptr));
      assert(!ev->write_other_ptr->write_other_ptr);
      ev->write_other_ptr->write_other_ptr = ev;
      assert(ev->write_other_ptr->trace_id() == ev->write_other_trace_id);
      assert(ev->write_other_ptr->write_other_trace_id == ev->trace_id());
      store_buffer[ev->thread_id()][ev->ml].pop_front();
      // WB -> WM thread order
      assert(!po.has_edge(ev, ev->write_other_ptr));
      if (!po.has_edge(ev->write_other_ptr, ev))
        po.add_edge(ev->write_other_ptr, ev);
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
    //
    if (is_read(ev)) {
      // Check if nd is not annotated yet
      if (!annotationPtr || !(annotationPtr->defines(ev))) {
        // This will be the last node for the corresponding thread
        has_unannotated_read_or_lock.insert(ev->thread_id());
      } else {
        assert(annotationPtr->defines(ev));
      }
      // Cache - wm
      if (!cache.wm.count(ev->ml))
        cache.wm.emplace
          (ev->ml, std::unordered_map<unsigned, std::vector<const ZEvent *>>());
      // Cache - readWB
      assert(!cache.readWB.count(ev));
      if (last_bwrite.count(ev->thread_id()) &&
          last_bwrite[ev->thread_id()].count(ev->ml)) {
        assert(is_writeB(last_bwrite[ev->thread_id()][ev->ml]));
        cache.readWB.emplace(ev, last_bwrite[ev->thread_id()][ev->ml]);
      }
      else
        cache.readWB.emplace(ev, nullptr);
    }

    // Handle Mutex events
    //
    if (is_mutex_init(ev)) {
      assert(!mutex_inits.count(ev->ml));
      mutex_inits.emplace(ev->ml, ev);
    }
    if (is_mutex_destroy(ev)) {
      assert(!mutex_destroys.count(ev->ml));
      mutex_destroys.emplace(ev->ml, ev);
    }
    if (is_lock(ev)) {
      // Check if ev is not annotated yet
      if (!annotationPtr || !(annotationPtr->location_has_some_lock(ev))) {
        // No annotated lock has held this location yet
        // This will be the last node for the corresponding thread
        has_unannotated_read_or_lock.insert(ev->thread_id());
      } else {
        // Some lock has already happened on this location
        if (annotationPtr->is_last_lock(ev)) {
          // This is the last annotated lock for this location
          assert(!found_last_lock_for_location.count(ev->ml));
          found_last_lock_for_location.insert(ev->ml);
        } else if (found_last_lock_for_location.count(ev->ml)) {
          // This is a lock after the last annotated lock for this location
          // This will be the last node for the corresponding thread
          has_unannotated_read_or_lock.insert(ev->thread_id());
        }
      }
    }
    if (is_lock(ev) || is_unlock(ev)) {
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
    auto thr_added = getThreadID(spwn->childs_cpid.get_proc_seq());
    assert(!thr_added.second);
    for (int aux : auxes(thr_added.first)) {
      assert(!(*this)(thr_added.first, aux).empty());
      const ZEvent *nthr = event(thr_added.first, aux, 0);

      assert(!po.has_edge(nthr, spwn));
      if (!po.has_edge(spwn, nthr))
        po.add_edge(spwn, nthr);
    }
  }

  // EDGES - joins
  for (const ZEvent *jn : joins) {
    auto thr_joined = getThreadID(jn->childs_cpid.get_proc_seq());
    assert(!thr_joined.second);
    for (int aux : auxes(thr_joined.first)) {
      assert(!(*this)(thr_joined.first, aux).empty());
      unsigned lastev_idx = (*this)(thr_joined.first, aux).size() - 1;
      const ZEvent *wthr = event(thr_joined.first, aux, lastev_idx);

      assert(!po.has_edge(jn, wthr));
      if (!po.has_edge(wthr, jn))
        po.add_edge(wthr, jn);
    }
  }

  // EDGES - mutex inits
  for (auto& in : mutex_inits) {
    const SymAddrSize& loc_init = in.first;
    const ZEvent *ev_init = in.second;

    for (auto& tid_ev_first : mutex_first[loc_init]) {
      const ZEvent *ev_first = tid_ev_first.second;

      assert(!po.has_edge(ev_first, ev_init));
      if (!po.has_edge(ev_init, ev_first))
        po.add_edge(ev_init, ev_first);
    }
  }

  // EDGES - mutex destroys
  for (auto& de : mutex_destroys) {
    const SymAddrSize& loc_destroy = de.first;
    const ZEvent *ev_destroy = de.second;

    for (auto& tid_ev_last : mutex_last[loc_destroy]) {
      const ZEvent *ev_last = tid_ev_last.second;

      assert(!po.has_edge(ev_destroy, ev_last));
      if (!po.has_edge(ev_last, ev_destroy))
        po.add_edge(ev_last, ev_destroy);
    }
  }
}


/* *************************** */
/* MAIN ALGORITHM              */
/* *************************** */

int ZGraph::getTailWindex(const SymAddrSize& ml, unsigned thr, int evX) const
{
  assert(thr < number_of_threads());
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


const ZEvent * ZGraph::getTailW
(const SymAddrSize& ml, unsigned thr, int evX) const
{
  int idx = getTailWindex(ml, thr, evX);
  if (idx == -1)
    return nullptr;

  assert(idx < (int) cache.wm.at(ml).at(thr).size());
  auto res = cache.wm.at(ml).at(thr)[idx];
  assert(is_writeM(res) && res->ml == ml &&
         res->thread_id() == thr && res->aux_id() != -1 &&
         evX >= 0 && res->event_id() <= (unsigned) evX);
  return res;
}


int ZGraph::getLatestNotAfterIndex(const ZEvent *read, unsigned thr, const ZPartialOrder& partial) const
{
  assert(is_read(read));
  assert(has_event(read));
  assert(thr < number_of_threads());
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


const ZEvent * ZGraph::getLatestNotAfter(const ZEvent *read, unsigned thr, const ZPartialOrder& partial) const
{
  assert(is_read(read));
  assert(has_event(read));
  int idx = getLatestNotAfterIndex(read, thr, partial);
  if (idx == -1)
    return nullptr;

  assert(idx < (int) cache.wm.at(read->ml).at(thr).size());
  auto res = cache.wm.at(read->ml).at(thr)[idx];
  assert(is_writeM(res) && same_ml(res, read) &&
         res->thread_id() == thr && res->aux_id() != -1 &&
         !partial.has_edge(read, res));
  return res;
}


const ZEvent * ZGraph::getLocalBufferW(const ZEvent *read) const
{
  assert(is_read(read));
  assert(has_event(read));
  assert(cache.readWB.count(read));
  const ZEvent *local = cache.readWB.at(read);
  assert(!local || is_writeB(local));
  return local;
}


std::list<const ZEvent *> ZGraph::eventsToMutate(const ZAnnotation& annotation) const
{
  auto res = std::list<const ZEvent *>();
  unsigned thr_id = 0;
  while (hasThreadAux(thr_id, -1)) {
    assert(!(*this)(thr_id, -1).empty());
    auto lastEv = event(thr_id, -1, (*this)(thr_id, -1).size() - 1);
    if ((is_read(lastEv) && !annotation.defines(lastEv)) ||
        (is_lock(lastEv) && !annotation.is_last_lock(lastEv)))
      res.push_back(lastEv);
    thr_id++;
  }

  // Sort the events based on a specified order
  res.sort(ZEventPtrComp());

  return res;
}


std::list<ZObs> ZGraph::getObsCandidates
(const ZEvent *read, const ZAnnotationNeg& negative) const
{
  assert(is_read(read));
  assert(has_event(read));

  std::list<ZObs> res;

  std::unordered_set<const ZEvent *> notCovered;
  std::unordered_set<const ZEvent *> mayBeCovered;
  // From the value (evid) and below, everything is covered from read by some other write
  std::unordered_map
    <std::pair<unsigned, int>, int> covered;
  for (unsigned covthr = 0; covthr < number_of_threads(); ++covthr)
    for (int covaux : auxes(covthr))
      covered.emplace(std::pair<unsigned, int>(covthr, covaux), -1);

  // Handle other threads
  for (unsigned tid = 0; tid < number_of_threads(); ++tid) {
    int auxid = auxForMl(read->ml, tid);
    if (read->thread_id() != tid && auxid >= 0) {
      int su = po.succ(read, tid, auxid).second;
      if (su >= (int) (*this)(tid, auxid).size()) {
        assert(su == INT_MAX);
        su = (int) (*this)(tid, auxid).size();
      }
      // Dumb
      /*
      int curr = su - 1;
      while (curr >= 0) {
        const ZEvent *rem = (*this)(tid, auxid)[curr];
        assert(is_writeM(rem) && !po.has_edge(read, rem));
        if (po.has_edge(rem, read)) {
          if (same_ml(read, rem)) {
            mayBeCovered.emplace(rem);
            // Update cover
            for (unsigned covthr = 0; covthr < number_of_threads(); ++covthr) {
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
          if (same_ml(read, rem))
            notCovered.emplace(rem);
          curr--;
        }
      }
      */
      // TailW index to cache.wm
      int wm_idx = getTailWindex(read->ml, tid, su - 1);
      assert(wm_idx == getLatestNotAfterIndex(read, tid, po));
      while (wm_idx > -1) {
        const ZEvent *rem = cache.wm.at(read->ml).at(tid).at(wm_idx);
        assert(rem && is_writeM(rem) && same_ml(read, rem));
        assert(!po.has_edge(read, rem));
        if (po.has_edge(rem, read)) {
          // Others in this thread are covered from read by rem
          assert((int) rem->event_id() <= po.pred(read, tid, auxid).second);
          mayBeCovered.emplace(rem);
          // Update cover
          for (unsigned covthr = 0; covthr < number_of_threads(); ++covthr) {
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

  // Handle thread of read
  const ZEvent *localB = cache.readWB.at(read);
  if (localB) {
    // There is a local write for read
    assert(is_writeB(localB) && same_ml(localB, read) &&
           localB->thread_id() == read->thread_id() &&
           localB->aux_id() == -1 &&
           localB->event_id() < read->event_id());
    const ZEvent *localM = localB->write_other_ptr;
    // Update cover caused by localB
    for (unsigned covthr = 0; covthr < number_of_threads(); ++covthr) {
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
    assert(is_writeM(localM) && same_ml(localB, localM) &&
           po.has_edge(localB, localM));
    for (auto it = notCovered.begin(); it != notCovered.end(); ) {
      const ZEvent *rem = *it;
      assert(is_writeM(rem) && !po.are_ordered(read, rem));
      if (po.has_edge(rem, localM))
        it = notCovered.erase(it);
      else
        ++it;
    }
    for (auto it = mayBeCovered.begin(); it != mayBeCovered.end(); ) {
      const ZEvent *rem = *it;
      assert(is_writeM(rem) && po.has_edge(rem, read));
      if (po.has_edge(rem, localM))
        it = mayBeCovered.erase(it);
      else
        ++it;
    }
    // Check if not covered by some remote
    bool localCovered = false;
    for (const auto& rem : mayBeCovered) {
      assert(is_writeM(rem) && po.has_edge(rem, read));
      if (po.has_edge(localM, rem)) {
        localCovered = true;
        break;
      }
    }
    // Add obs if not covered and not forbidden
    if (!localCovered && !negative.forbids(read, localB))
      res.emplace_back(localB);
  } else {
    // No local write for read
    if (mayBeCovered.empty()) {
      // Consider initial event if not forbidden
      if (!negative.forbids_initial(read))
        res.emplace_back(INT_MAX, INT_MAX);
    }
  }

  // Take candidates unordered with read
  for (const auto& remM : notCovered) {
    assert(is_writeM(remM));
    const ZEvent *remB = remM->write_other_ptr;
    assert(is_writeB(remB) && same_ml(remB, remM) &&
           po.has_edge(remB, remM));
    // Add if not forbidden
    if (!negative.forbids(read, remB))
      res.emplace_back(remB);
  }

  // Take candidates that happen before read
  for (const auto& remM : mayBeCovered) {
    assert(is_writeM(remM));
    auto covta = std::pair<unsigned, int>
      (remM->thread_id(), remM->aux_id());
    assert(covered.count(covta));
    if ((int) remM->event_id() > covered[covta]) {
      // Not covered
      const ZEvent *remB = remM->write_other_ptr;
      assert(is_writeB(remB) && same_ml(remB, remM) &&
             po.has_edge(remB, remM));
      // Add if not forbidden
      if (!negative.forbids(read, remB))
        res.emplace_back(remB);
    }
  }

  res.sort();
  return res;
}
