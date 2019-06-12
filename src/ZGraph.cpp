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

#include <list>

#include "ZHelpers.h"
#include "ZGraph.h"
#include "ZBuilderTSO.h"


// Empty
ZGraph::ZGraph()
  : basis(*this), po(this->basis), cache()
{
  assert(&(basis.graph) == this);
  assert(&(po.basis) == &basis);
  assert(empty());
}


// Initial
ZGraph::ZGraph(const std::vector<ZEvent>& trace, int star_root_index)
  : basis(*this, star_root_index),
    po(this->basis),
    cache()
{
  assert(&(basis.graph) == this);
  assert(&(po.basis) == &basis);
  traceToPO(trace, nullptr);
  assert(basis.root() < basis.size() &&
         "Root index too big (not enough threads in the initial trace)");
}


// Moving
ZGraph::ZGraph(ZGraph&& oth)
  : basis(std::move(oth.basis)),
    po(std::move(oth.po)),
    cache(std::move(oth.cache))
{
  assert(oth.empty());
}


// Partial order that will be moved as original
// Trace and annotation that will extend this copy of the graph
ZGraph::ZGraph(const ZGraph& oth,
       ZPartialOrder&& po,
       const std::vector<ZEvent>& trace,
       const ZAnnotation& annotation)
  : basis(oth.basis, *this),
    po(po, this->basis),
    cache()
{
  assert(&(basis.graph) == this);
  assert(&(po.basis) == &basis);
  traceToPO(trace, &annotation);
}


/* *************************** */
/* GRAPH EXTENSION             */
/* *************************** */

// Extends this graph so it corresponds to 'trace'
// Check the header file for the method description
void ZGraph::traceToPO(const std::vector<ZEvent>& trace,
                       const ZAnnotation *annotationPtr)
{
  assert(!trace.empty());
  assert(trace.size() >= basis.events_size());

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
  proc_seq_within_po.insert(trace[0].cpid.get_proc_seq());

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
    if (!proc_seq_within_po.count(ev->cpid.get_proc_seq())) {
      // This is a forbidden process because it
      // is created after a so-far unannotated node
      forbidden_processes_ipid.insert(ev->iid.get_pid());
      continue;
    }

    unsigned thr_idx = INT_MAX;

    // Check if this process is already known
    auto ipidit = ipid_to_thraux.find(ev->iid.get_pid());
    if (ipidit != ipid_to_thraux.end()) {
      thr_idx = ipidit->second.first;
      assert(ev->auxID() == ipidit->second.second);
    } else {
      auto thr_new = basis.getThreadID(ev);
      thr_idx = thr_new.first;
      // add to ipid cache for faster lookup next time
      ipid_to_thraux.emplace(ev->iid.get_pid(),
                             std::pair<unsigned, int>(thr_idx,
                                                      ev->auxID()));
    }

    ev->_thread_id = thr_idx;
    auto thraux = std::pair<unsigned, int>(ev->threadID(), ev->auxID());

    // Check possible cases why we cannot add this event
    //
    // Check if we already haven't added a thraux-predecessor of this event
    if (thraux_with_event_we_dont_add.count(thraux)) {
      assert(!basis.hasThreadAux(thraux) ||
             cur_evidx.at(thraux) == basis(thraux).size());
      continue;
    }
    // Check if thraux already has an unannotated read
    if (has_unannotated_read_or_lock.count(thraux.first) &&
        thraux.second == -1) {
      // This event happens after something we need to
      // annotate first, so we don't add it into the graph
      assert(!basis.hasThreadAux(thraux) ||
             cur_evidx.at(thraux) == basis(thraux).size());
      thraux_with_event_we_dont_add.insert(thraux);
      continue;
    }
    // Joins of threads with an event we didn't include also
    // can not be included and neither any of their successors
    if (isJoin(ev)) {
      bool dontInclude = false;
      auto childs_proc_it = proc_seq_within_po.find(ev->childs_cpid.get_proc_seq());
      if (childs_proc_it == proc_seq_within_po.end()) {
        dontInclude = true;
      } else {
        auto thr_added = basis.getThreadID(ev->childs_cpid.get_proc_seq());
        assert(!thr_added.second);
        auto childs_thraux = std::pair<unsigned, int>(thr_added.first, -1);
        // It is enough to check the presence of the real thread
        if (thraux_with_event_we_dont_add.count(childs_thraux))
          dontInclude = true;
      }
      if (dontInclude) {
        assert(!basis.hasThreadAux(thraux) ||
               cur_evidx.at(thraux) == basis(thraux).size());
        thraux_with_event_we_dont_add.insert(thraux);
        continue;
      }
    }
    // Memory-writes cannot be added if their buffer-write
    // was not added
    if (isWriteM(ev)) {
      assert(ev->auxID() != -1 && "Only auxiliary threads");
      if (!store_buffer.count(ev->threadID()) ||
          !store_buffer[ev->threadID()].count(ev->ml) ||
          store_buffer[ev->threadID()][ev->ml].empty()) {
        assert(!basis.hasThreadAux(thraux) ||
               cur_evidx.at(thraux) == basis(thraux).size());
        thraux_with_event_we_dont_add.insert(thraux);
        continue;
      }
    }

    // We will add the event
    // First, add line if not already present
    if (!basis.hasThreadAux(thraux)) {
      // A new thread/aux appeared,
      // we've checked above that it is allowed
      assert(!cur_evidx.count(thraux));
      cur_evidx.emplace(thraux, 0);

      basis.addLine(ev);
      po.addLine(ev);
    }

    unsigned ev_idx = cur_evidx[thraux];

    assert(ev->eventID() == ev_idx);
    if (ev_idx == basis(thraux).size()) {
      // New event
      basis.addEvent(ev);
      po.addEvent(ev);

      // XXX: function pointer calls not handled
      if (isSpawn(ev))
        spawns.push_back(ev);
      if (isJoin(ev))
        joins.push_back(ev);

    } else {
      // Already known event
      const ZEvent *oldev = basis.getEvent(ev->threadID(),
                                           ev->auxID(),
                                           ev->eventID());
      basis.replaceEvent(oldev, ev);
    }

    ++cur_evidx[thraux];

    // Handle Fence
    if (ev->fence) {
      if (store_buffer.count(ev->threadID()))
        for (const auto& ml_last : last_mwrite[ev->threadID()]) {
          auto last = ml_last.second;
          assert(isWriteM(last));
          assert(!po.hasEdge(ev, last));
          if (!po.hasEdge(last, ev))
            po.addEdge(last, ev);
        }
    }
    // Handle Spawn
    if (isSpawn(ev)) {
      proc_seq_within_po.insert(ev->childs_cpid.get_proc_seq());
      assert(basis.number_of_threads() > 0);
      if (ev->fence && basis.number_of_threads() == 1) {
        // Only one thread now, and everything there
        // happens before this spawn event, so nothing
        // in this thread has to be chrono ordered
        cache.chrono.clear();
      }
    }
    // Handle Buffer-Write
    if (isWriteB(ev)) {
      assert(ev->auxID() == -1 && "Only real threads");
      // Store buffer
      if (!store_buffer.count(ev->threadID()))
        store_buffer.emplace
          (ev->threadID(),
           std::unordered_map
           <SymAddrSize, std::list<const ZEvent *>>());
      if (!store_buffer[ev->threadID()].count(ev->ml))
        store_buffer[ev->threadID()].emplace(ev->ml,
                                             std::list<const ZEvent *>());
      store_buffer[ev->threadID()][ev->ml].push_back(ev);
      // Last Bwrite
      if (!last_bwrite.count(ev->threadID()))
        last_bwrite.emplace
          (ev->threadID(), std::unordered_map<SymAddrSize, const ZEvent *>());
      auto it = last_bwrite[ev->threadID()].find(ev->ml);
      if (it != last_bwrite[ev->threadID()].end())
        it = last_bwrite[ev->threadID()].erase(it);
      last_bwrite[ev->threadID()].emplace_hint
        (it, ev->ml, ev);
    }
    // Handle Memory-Write
    if (isWriteM(ev)) {
      assert(ev->auxID() != -1 && "Only auxiliary threads");
      // Store buffer
      assert(store_buffer.count(ev->threadID()) &&
             store_buffer[ev->threadID()].count(ev->ml));
      assert(!ev->write_other_ptr);
      ev->write_other_ptr = store_buffer[ev->threadID()][ev->ml].front();
      assert(isWriteB(ev->write_other_ptr) && sameMl(ev, ev->write_other_ptr));
      assert(!ev->write_other_ptr->write_other_ptr);
      ev->write_other_ptr->write_other_ptr = ev;
      assert(ev->write_other_ptr->traceID() == ev->write_other_trace_id);
      assert(ev->write_other_ptr->write_other_trace_id == ev->traceID());
      store_buffer[ev->threadID()][ev->ml].pop_front();
      // WB -> WM thread order
      assert(!po.hasEdge(ev, ev->write_other_ptr));
      if (!po.hasEdge(ev->write_other_ptr, ev))
        po.addEdge(ev->write_other_ptr, ev);
      // Last Mwrite
      if (!last_mwrite.count(ev->threadID()))
        last_mwrite.emplace
          (ev->threadID(), std::unordered_map<SymAddrSize, const ZEvent *>());
      auto it = last_mwrite[ev->threadID()].find(ev->ml);
      if (it != last_mwrite[ev->threadID()].end())
        it = last_mwrite[ev->threadID()].erase(it);
      last_mwrite[ev->threadID()].emplace_hint
        (it, ev->ml, ev);
      // Cache - wm
      if (!cache.wm.count(ev->ml))
        cache.wm.emplace
          (ev->ml, std::unordered_map<unsigned, std::vector<const ZEvent *>>());
      if (!cache.wm[ev->ml].count(ev->threadID()))
        cache.wm[ev->ml].emplace(ev->threadID(), std::vector<const ZEvent *>());
      cache.wm[ev->ml][ev->threadID()].push_back(ev);
      // Cache - chrono
      if (!basis.isRoot(ev)) {
        if (!cache.chrono.count(ev->threadID()))
          cache.chrono.emplace
            (ev->threadID(), std::list<const ZEvent *>());
        cache.chrono[ev->threadID()].push_back(ev);
      }
    }
    // Handle Read
    if (isRead(ev)) {
      // Check if nd is not annotated yet
      if (!annotationPtr || !(annotationPtr->defines(ev))) {
        // This will be the last node for the corresponding thread
        has_unannotated_read_or_lock.insert(ev->threadID());
      }
      // Cache - wm
      if (!cache.wm.count(ev->ml))
        cache.wm.emplace
          (ev->ml, std::unordered_map<unsigned, std::vector<const ZEvent *>>());
      // Cache - readWB
      assert(!cache.readWB.count(ev));
      if (last_bwrite.count(ev->threadID()) &&
          last_bwrite[ev->threadID()].count(ev->ml)) {
        assert(isWriteB(last_bwrite[ev->threadID()][ev->ml]));
        cache.readWB.emplace(ev, last_bwrite[ev->threadID()][ev->ml]);
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
    if (isLock(ev)) {
      // Check if ev is not annotated yet
      if (!annotationPtr || !(annotationPtr->locationHasSomeLock(ev))) {
        // No annotated lock has held this location yet
        // This will be the last node for the corresponding thread
        has_unannotated_read_or_lock.insert(ev->threadID());
      } else {
        // Some lock has already happened on this location
        if (annotationPtr->isLastLock(ev)) {
          // This is the last annotated lock for this location
          assert(!found_last_lock_for_location.count(ev->ml));
          found_last_lock_for_location.insert(ev->ml);
        } else if (found_last_lock_for_location.count(ev->ml)) {
          // This is a lock after the last annotated lock for this location
          // This will be the last node for the corresponding thread
          has_unannotated_read_or_lock.insert(ev->threadID());
        }
      }
    }
    if (isLock(ev) || isUnlock(ev)) {
      // For each ml, for each thread, remember first access
      if (!mutex_first.count(ev->ml)) {
        mutex_first.emplace(ev->ml, std::unordered_map<unsigned, const ZEvent *>());
        mutex_first[ev->ml].reserve(8);
      }
      if (!mutex_first[ev->ml].count(ev->threadID()))
        mutex_first[ev->ml].emplace(ev->threadID(), ev);
      // For each ml, for each thread, remember last access
      if (!mutex_last.count(ev->ml)) {
        mutex_last.emplace(ev->ml, std::unordered_map<unsigned, const ZEvent *>());
        mutex_last[ev->ml].reserve(8);
      }
      mutex_last[ev->ml][ev->threadID()] = ev;
    }
  } // end of loop for traversing trace and creating nodes

  basis.shrink();
  po.shrink();

  #ifndef NDEBUG
  assert(basis.size() == cur_evidx.size());
  for (auto& thraux_numevents : cur_evidx) {
    assert(thraux_numevents.second == basis(thraux_numevents.first).size()
           && "Didn't go through entire original part of the basis");
  }
  #endif

  // EDGES - spawns
  for (const ZEvent *spwn : spawns) {
    auto thr_added = basis.getThreadID(spwn->childs_cpid.get_proc_seq());
    assert(!thr_added.second);
    int aux = -1;
    while (basis.hasThreadAux(thr_added.first, aux)) {
      assert(!basis(thr_added.first, aux).empty());
      const ZEvent *nthr = basis.getEvent(thr_added.first, aux, 0);

      assert(!po.hasEdge(nthr, spwn));
      if (!po.hasEdge(spwn, nthr))
        po.addEdge(spwn, nthr);
      aux++;
    }
  }

  // EDGES - joins
  for (const ZEvent *jn : joins) {
    auto thr_joined = basis.getThreadID(jn->childs_cpid.get_proc_seq());
    assert(!thr_joined.second);
    int aux = -1;
    while (basis.hasThreadAux(thr_joined.first, aux)) {
      assert(!basis(thr_joined.first, aux).empty());
      unsigned lastev_idx = basis(thr_joined.first, aux).size() - 1;
      const ZEvent *wthr = basis.getEvent(thr_joined.first, aux, lastev_idx);

      assert(!po.hasEdge(jn, wthr));
      if (!po.hasEdge(wthr, jn))
        po.addEdge(wthr, jn);

      aux++;
    }
    assert(aux >= 0);
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
}


/* *************************** */
/* MAIN ALGORITHM              */
/* *************************** */

int ZGraph::getTailWindex(const SymAddrSize& ml, unsigned thr, unsigned ev) const
{
  assert(thr < basis.number_of_threads());
  assert(cache.wm.count(ml));
  if (!cache.wm.at(ml).count(thr))
    return -1;

  const auto& writes = cache.wm.at(ml).at(thr);
  if (writes.empty())
    return -1;

  int low = 0;
  if (ev < writes[low]->eventID())
    return -1;

  int high = writes.size() - 1;
  if (ev >= writes[high]->eventID())
    return high;

  assert(low < high);
  if (low + 1 == high) {
    assert(ev >= writes[low]->eventID() &&
           ev < writes[high]->eventID());
    return low;
  }

  // Low represents something that
  // can possibly be the answer
  // High represents something that
  // is above the answer
  // Do binary search
  while (true) {
    assert(low + 1 < high);
    assert(ev >= writes[low]->eventID() &&
           ev < writes[high]->eventID());
    int mid = ((high - low) / 2) + low;
    assert(low < mid && mid < high);

    if (ev >= writes[mid]->eventID())
      low = mid;
    else
      high = mid;

    if (low + 1 == high) {
      assert(ev >= writes[low]->eventID() &&
             ev < writes[high]->eventID());
      return low;
    }
  }
}


const ZEvent * ZGraph::getTailW
(const SymAddrSize& ml, unsigned thr, unsigned ev) const
{
  int idx = getTailWindex(ml, thr, ev);
  if (idx == -1)
    return nullptr;

  assert(idx < (int) cache.wm.at(ml).at(thr).size());
  auto res = cache.wm.at(ml).at(thr)[idx];
  assert(isWriteM(res) && res->ml == ml &&
         res->threadID() == thr && res->auxID() != -1 &&
         res->eventID() <= ev);
  return res;
}


int ZGraph::getLatestNotAfterIndex(const ZEvent *read, unsigned thr, const ZPartialOrder& partial) const
{
  assert(isRead(read));
  assert(basis.hasEvent(read));
  assert(thr < basis.number_of_threads());
  assert(thr != read->threadID());
  assert(cache.wm.count(read->ml));

  int aux = basis.auxForMl(read->ml, thr);
  if (aux == -1) {
    assert(!cache.wm.at(read->ml).count(thr));
    return -1;
  }
  assert(aux == 0 || cache.wm.at(read->ml).count(thr));

  int su = partial.succ(read, thr, aux).second;
  if (su >= (int) basis(thr, aux).size()) {
    assert(su == INT_MAX);
    su = (int) basis(thr, aux).size();
  }
  // TailW index to cache.wm
  return getTailWindex(read->ml, thr, su - 1);
}


const ZEvent * ZGraph::getLatestNotAfter(const ZEvent *read, unsigned thr, const ZPartialOrder& partial) const
{
  int idx = getLatestNotAfterIndex(read, thr, partial);
  if (idx == -1)
    return nullptr;

  assert(idx < (int) cache.wm.at(read->ml).at(thr).size());
  auto res = cache.wm.at(read->ml).at(thr)[idx];
  assert(isWriteM(res) && sameMl(res, read) &&
         res->threadID() == thr && res->auxID() != -1 &&
         !partial.hasEdge(read, res));
  return res;
}


const ZEvent * ZGraph::getLocalBufferW(const ZEvent *read) const
{
  assert(isRead(read));
  assert(basis.hasEvent(read));
  assert(cache.readWB.count(read));
  const ZEvent *local = cache.readWB.at(read);
  assert(!local || isWriteB(local));

  return local;
}


std::list<const ZEvent *> ZGraph::getEventsToMutate(const ZAnnotation& annotation) const
{
  auto res = std::list<const ZEvent *>();
  unsigned thr_id = 0;
  while (basis.hasThreadAux(thr_id, -1)) {
    assert(!basis(thr_id, -1).empty());
    auto lastEv = basis.getEvent(thr_id, -1,
                                 basis(thr_id, -1).size() - 1);
    if ((isRead(lastEv) && !annotation.defines(lastEv)) ||
        (isLock(lastEv) && !annotation.isLastLock(lastEv)))
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
  assert(isRead(read));
  assert(basis.hasEvent(read));

  std::list<ZObs> res;

  std::unordered_set<const ZEvent *> notCovered;
  std::unordered_set<const ZEvent *> mayBeCovered;
  // From the value (evid) and below, everything is covered from read by some other write
  std::unordered_map
    <std::pair<unsigned, int>, int> covered;
  for (unsigned covthr = 0; covthr < basis.number_of_threads(); ++covthr)
    for (int covaux : basis.auxes(covthr))
      covered.emplace(std::pair<unsigned, int>(covthr, covaux), -1);

  // Handle other threads
  for (unsigned tid = 0; tid < basis.number_of_threads(); ++tid) {
    int auxid = basis.auxForMl(read->ml, tid);
    if (read->threadID() != tid && auxid >= 0) {
      int su = po.succ(read, tid, auxid).second;
      if (su >= (int) basis(tid, auxid).size()) {
        assert(su == INT_MAX);
        su = (int) basis(tid, auxid).size();
      }
      // TailW index to cache.wm
      int wm_idx = getTailWindex(read->ml, tid, su - 1);
      assert(wm_idx == getLatestNotAfterIndex(read, tid, po));
      while (wm_idx > -1) {
        const ZEvent *rem = cache.wm.at(read->ml).at(tid).at(wm_idx);
        assert(rem && isWriteM(rem) && sameMl(read, rem));
        assert(!po.hasEdge(read, rem));
        if (po.hasEdge(rem, read)) {
          // Others in this thread are covered from read by rem
          assert((int) rem->eventID() <= po.pred(read, tid, auxid).second);
          mayBeCovered.emplace(rem);
          // Update cover
          for (unsigned covthr = 0; covthr < basis.number_of_threads(); ++covthr) {
            if (covthr != rem->threadID()) {
              for (int covaux : basis.auxes(covthr)) {
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
    assert(isWriteB(localB) && sameMl(localB, read) &&
           localB->threadID() == read->threadID() &&
           localB->auxID() == -1 &&
           localB->eventID() < read->eventID());
    // Update cover caused by localB
    for (unsigned covthr = 0; covthr < basis.number_of_threads(); ++covthr) {
      if (covthr != localB->threadID()) {
        for (int covaux : basis.auxes(covthr)) {
          int newcov = po.pred(localB, covthr, covaux).second;
          auto covta = std::pair<unsigned, int>(covthr, covaux);
          assert(covered.count(covta));
          if (newcov > covered[covta])
            covered[covta] = newcov;
        }
      }
    }
    // Delete remotes that happen after localM
    const ZEvent *localM = localB->write_other_ptr;
    assert(isWriteM(localM) && sameMl(localB, localM) &&
           po.hasEdge(localB, localM));
    for (auto it = notCovered.begin(); it != notCovered.end(); ) {
      const ZEvent *rem = *it;
      assert(isWriteM(rem) && !po.areOrdered(read, rem));
      if (po.hasEdge(localM, rem))
        it = notCovered.erase(it);
      else
        ++it;
    }
    for (auto it = mayBeCovered.begin(); it != mayBeCovered.end(); ) {
      const ZEvent *rem = *it;
      assert(isWriteM(rem) && po.hasEdge(rem, read));
      if (po.hasEdge(localM, rem))
        it = mayBeCovered.erase(it);
      else
        ++it;
    }
    // Add obs if not forbidden
    if (!negative.forbids(read, localB))
      res.emplace_back(localB);
  } else {
    // No local write for read
    if (mayBeCovered.empty()) {
      // Consider initial event if not forbidden
      if (!negative.forbidsInitialEvent(read))
        res.emplace_back(INT_MAX, INT_MAX);
    }
  }

  // Take candidates unordered with read
  for (const auto& remM : notCovered) {
    assert(isWriteM(remM));
    const ZEvent *remB = remM->write_other_ptr;
    assert(isWriteB(remB) && sameMl(remB, remM) &&
           po.hasEdge(remB, remM));
    // Add if not forbidden
    if (!negative.forbids(read, remB))
      res.emplace_back(remB);
  }

  // Take candidates that happen before read
  for (const auto& remM : mayBeCovered) {
    assert(isWriteM(remM));
    auto covta = std::pair<unsigned, int>
      (remM->threadID(), remM->auxID());
    assert(covered.count(covta));
    if ((int) remM->eventID() > covered[covta]) {
      // Not covered
      const ZEvent *remB = remM->write_other_ptr;
      assert(isWriteB(remB) && sameMl(remB, remM) &&
             po.hasEdge(remB, remM));
      // Add if not forbidden
      if (!negative.forbids(read, remB))
        res.emplace_back(remB);
      }
  }

  res.sort();
  return res;
}


/*
const Node * ZGraph::getTailWcandidate(const Node *nd, unsigned thr_id,
                                              const PartialOrder& po) const
{
  assert(isRead(nd));
  const ThreadPairsVclocks& succ = *(po.first);

  // Get index where to start the search from
  int ev_id = (nd->getProcessID() == thr_id) ? nd->getEventID() - 1
            : succ[nd->getProcessID()][thr_id][nd->getEventID()] - 1;
  if (ev_id == INT_MAX - 1)
    ev_id = processes[thr_id].size() - 1;
  assert(ev_id < (int) processes[thr_id].size());

  if (ev_id == -1) // There are no writes in thr_id not happening after nd
    return (nd->getProcessID() == thr_id)
      ? initial_node : nullptr; // Initial node treated as from same thread

  // TAIL WRITE CANDIDATES CACHE
  // [ml][tid][evid] returns idx of first event of thread-tid writing to ml
  // starting from AND INCLUDING evid and going back - (evid, evid-1, .., 0)
  // returns -1 if there is no such write

  auto ml_cache = tw_candidate.find(nd->getEvent()->ml);
  assert(ml_cache != tw_candidate.end()
         && "Cache not set up for the ml of this read");
  assert(thr_id < ml_cache->second.size() && ev_id >= 0 &&
         ev_id < (int) ml_cache->second[thr_id].size());

  // Use tail write candidates cache to get the result
  int tw_evidx = ml_cache->second[thr_id][ev_id];
  if (tw_evidx == -1)
    return (nd->getProcessID() == thr_id)
      ? initial_node : nullptr; // Initial node treated as from same thread
  assert(tw_evidx >= 0 && tw_evidx <= ev_id);
  const Node *result = processes[thr_id][tw_evidx];
  assert(isWrite(result) && sameMl(result, nd));
  return result;
}

std::pair<const Node *, std::unordered_set<const Node *>>
ZGraph::getTailWrites(const Node *nd, const PartialOrder& po) const
{
  assert(isRead(nd));

  // Get candidate for each thread
  auto result = std::unordered_set<const Node *>();
  for (unsigned thr_id = 0; thr_id < processes.size(); ++thr_id) {
    const Node *twcand = getTailWcandidate(nd, thr_id, po);
    if (twcand != nullptr)
      result.emplace(twcand);
  }

  // Retain those not happening before another candidate
  for (auto it = result.begin(); it != result.end(); ) {
    const Node *twcand = *it;
    bool tail = true;
    for (const Node *twother : result)
      if (twcand != twother && hasEdge(twcand, twother, po)) {
        tail = false; // HB another candidate, so not tail
        break;
      }
    if (!tail)
      it = result.erase(it);
    else
      ++it;
  }

  // Return root tail write separately
  const Node *rootTail = nullptr;
  for (auto it = result.begin(); it != result.end(); ) {
    if ((*it)->getProcessID() == starRoot() ||
        (nd->getProcessID() == starRoot() && *it == initial_node)) {
      assert(rootTail == nullptr);
      rootTail = *it;
      it = result.erase(it);
    } else
      ++it;
  }

  return {rootTail, result};
}

std::pair<const Node *, std::pair<int, int>>
ZGraph::getHeadWcandidate(const Node *nd, unsigned thr_id, const PartialOrder& po) const
{
  // TAIL WRITE CANDIDATES CACHE
  // [ml][tid][evid] returns idx of first event of thread-tid writing to ml
  // starting from AND INCLUDING evid and going back - (evid, evid-1, .., 0)
  // returns -1 if there is no such write

  if (nd->getProcessID() == thr_id) {
    // Looking for a HB candidate from the same thread
    if (nd->getEventID() == 0)
      return {initial_node, {0, -1}}; // Initial node treated as from same thread
    // Use cache to get the candidate happening before nd
    int ev_id = nd->getEventID() - 1;
    auto ml_cache = tw_candidate.find(nd->getEvent()->ml);
    assert(ml_cache != tw_candidate.end()
           && "Cache not set up for the ml of this read");
    assert(thr_id < ml_cache->second.size() && ev_id >= 0 &&
           ev_id < (int) ml_cache->second[thr_id].size());
    int tw_evidx = ml_cache->second[thr_id][ev_id];

    if (tw_evidx == -1)
      return {initial_node, {0, -1}}; // Initial node treated as from same thread
    assert(tw_evidx >= 0 && tw_evidx <= ev_id);
    const Node *before_nd = processes[thr_id][tw_evidx];
    assert(isWrite(before_nd) && sameMl(before_nd, nd));
    return {before_nd, {0, -1}};
  }

  // Looking for candidate from a different thread
  const ThreadPairsVclocks& succ = *(po.first);
  const ThreadPairsVclocks& pred = *(po.second);

  // Search for a head write candidate happening before nd
  const Node *before_nd = nullptr;
  int ev_id = pred[nd->getProcessID()][thr_id][nd->getEventID()];
  if (ev_id != -1) {
    assert (ev_id >= 0 && ev_id < (int) processes[thr_id].size());
    auto ml_cache = tw_candidate.find(nd->getEvent()->ml);
    assert(ml_cache != tw_candidate.end()
           && "Cache not set up for the ml of this read");
    assert(thr_id < ml_cache->second.size() && ev_id >= 0 &&
           ev_id < (int) ml_cache->second[thr_id].size());
    int tw_evidx = ml_cache->second[thr_id][ev_id];

    if (tw_evidx != -1) {
      // Got a head write candidate that happens before nd
      assert(tw_evidx >= 0 && tw_evidx <= ev_id);
      before_nd = processes[thr_id][tw_evidx];
      assert(isWrite(before_nd) && sameMl(before_nd, nd) &&
             hasEdge(before_nd, nd, po));
    }
  }
  // Return indices from-to for unordered head write candidate search
  int limit = succ[nd->getProcessID()][thr_id][nd->getEventID()];
  if (limit == INT_MAX)
    limit = processes[thr_id].size();
  return {before_nd, {ev_id+1, limit}};
}

std::pair<const Node *, std::unordered_set<const Node *>>
ZGraph::getHeadWrites(const Node *nd, const PartialOrder& po) const
{
  assert(isRead(nd));

  // Get before candidates and search indices for each thread
  auto befores = std::vector<const Node *>();
  befores.reserve(processes.size());
  auto search = std::vector<std::pair<int,int>>();
  search.reserve(processes.size());
  for (unsigned thr_id = 0; thr_id < processes.size(); ++thr_id) {
    auto before_search = getHeadWcandidate(nd, thr_id, po);
    befores.push_back(before_search.first);
    search.emplace_back(before_search.second.first,
                        before_search.second.second);
  }

  // Delete befores that are covered by another before
  for (unsigned thr_id = 0; thr_id < processes.size(); ++thr_id)
    if (befores[thr_id] != nullptr) {
      for (const Node *other_before : befores)
        if (other_before != nullptr && other_before != befores[thr_id] &&
            hasEdge(befores[thr_id], other_before, po)) {
          assert(hasEdge(other_before, nd, po));
          befores[thr_id] = nullptr;
          break;
        }
    }

  // All remaining befores are head writes
  // Threads without before can have a head write
  // that is unordered with nd, search for first
  // one like that (it can't be covered by any before)
  auto unords = std::vector<const Node *>(processes.size(), nullptr);
  for (unsigned thr_id = 0; thr_id < processes.size(); ++thr_id)
    if (befores[thr_id] == nullptr) {
      // Search for an unordered head candidate
      for (int unord_ev_id = search[thr_id].first;
           unord_ev_id < search[thr_id].second; ++unord_ev_id) {
        const Node *unord_nd = processes[thr_id][unord_ev_id];
        if (isWrite(unord_nd) && sameMl(unord_nd, nd)) {
          // Got a conflicting write unordered with nd, this is
          // the only possible unord head candidate for this thread
          assert(!areOrdered(nd, unord_nd, po));
          #ifndef NDEBUG
          // Assert it is not covered by any before
          // (because otherwise it would HB nd)
          for (const Node *before_nd : befores)
            assert(before_nd == nullptr ||
                   !hasEdge(unord_nd, before_nd, po));
          #endif
          unords[thr_id] = unord_nd;
          break;
        }
      }
    }

  // For every remaining unord: it is head write if
  // no remaining before or unord happens before it
  for (unsigned thr_id = 0; thr_id < processes.size(); ++thr_id)
    if (unords[thr_id] != nullptr) {
      bool head = true;
      for (const Node *other_before : befores)
        if (other_before != nullptr && hasEdge(other_before, unords[thr_id], po)) {
          // this unord is not head, because other_before happens before it
          head = false;
          break;
        }
      if (head)
        for (const Node *other_unord : unords)
          if (other_unord != nullptr && other_unord != unords[thr_id] &&
              hasEdge(other_unord, unords[thr_id], po)) {
            // this unord is not head, because other_unord happens before it
            head = false;
            break;
          }
      if (!head)
        unords[thr_id] = nullptr;
    }

  // Return root head write separately
  const Node *rootHead = nullptr;
  auto result = std::unordered_set<const Node *>();
  for (const Node *before : befores)
    if (before != nullptr) {
      if (before->getProcessID() == starRoot() ||
          (nd->getProcessID() == starRoot() && before == initial_node)) {
        assert(rootHead == nullptr);
        rootHead = before;
      } else
        result.emplace(before);
    }
  for (const Node *unord : unords)
    if (unord != nullptr) {
      if (unord->getProcessID() == starRoot() ||
          (nd->getProcessID() == starRoot() && unord == initial_node)) {
        assert(rootHead == nullptr);
        rootHead = unord;
      } else
        result.emplace(unord);
    }

  return {rootHead, result};
}

bool ZGraph::isObservable(const Node *nd, const PartialOrder& po) const
{
  assert(isWrite(nd));

  auto itmlR = readsRoot.find(nd->getEvent()->ml);
  if (itmlR != readsRoot.end())
    for (const Node *rootRead : itmlR->second) {
      bool observableBy = isObservableBy(nd, rootRead, po);
      if (observableBy)
        return true;
    }

  auto itmlN = readsNonroot.find(nd->getEvent()->ml);
  if (itmlN != readsNonroot.end())
    for (const Node *nonrootRead : itmlN->second) {
      bool observableBy = isObservableBy(nd, nonrootRead, po);
      if (observableBy)
        return true;
    }

  return false;
}

bool ZGraph::isObservableBy(const Node *writend, const Node *readnd,
                                   const PartialOrder& po) const
{
  assert(isWrite(writend) && isRead(readnd) &&
         sameMl(writend, readnd));

  if (hasEdge(readnd, writend, po)) {
    // Trivially unobservable
    return false;
  }

  if (!hasEdge(writend, readnd, po)) {
    // Trivially observable
    assert(!areOrdered(writend, readnd, po));
    return true;
  }

  assert(hasEdge(writend, readnd, po));
  // If a write happens before a read, it can be
  // observable by this read only if it is head
  auto heads = getHeadWrites(readnd, po);
  if (heads.first && heads.first == writend)
    return true;

  if (heads.second.size() > 0)
    for (const Node *nonrootHead : heads.second)
      if (nonrootHead == writend)
        return true;

  return false;
}

void ZGraph::orderEventMaz(const ZEvent *ev1, const ZAnnotation& annotation,
                                  bool newlyEverGoodWrite, const PartialOrder& po)
{
  assert(isWrite(ev1) || isRead(ev1));
  auto it = event_to_node.find(ev1);
  assert(it != event_to_node.end());
  const Node *nd1 = it->second;
  assert(nd1->getEvent() == ev1);
  assert(nd1->getProcessID() != starRoot());

  bool readOrEverGoodWrite = (isRead(ev1) || newlyEverGoodWrite ||
                              annotation.isEverGood(nd1));

  // Collect conflicting nonroot writes to order
  auto itnsw = wNonrootUnord.find(ev1->ml);
  assert(itnsw != wNonrootUnord.end());
  auto toOrder = std::vector<const Node *>();
  toOrder.reserve(itnsw->second.size());

  for (auto& writend : itnsw->second) {
    assert(writend->getProcessID() != starRoot());
    if (nd1->getProcessID() != writend->getProcessID()) {
      // If nd1 is a read, take all write candidates,
      // otherwise take only candidates such that
      // at least one of them is everGood
      if (readOrEverGoodWrite || annotation.isEverGood(writend))
        toOrder.push_back(writend);
    }
  }

  // If ev1 is write, also collect
  // conflicting annotated nonroot reads to order
  if (isWrite(ev1)) {
    auto itnsr = readsNonroot.find(ev1->ml);
    if (itnsr != readsNonroot.end()) {
      for (auto& readnd : itnsr->second)
        if (nd1->getProcessID() != readnd->getProcessID()
            && annotation.defines(readnd))
          toOrder.push_back(readnd);
    }
  }

  // Sort the vector so that the execution is deterministic
  if (toOrder.size() > 1) {
    auto comp = NodePtrComp();
    #ifndef NDEBUG
    for (auto it1 = toOrder.begin();
         it1 != toOrder.end(); ++it1)
      for (auto it2 = toOrder.begin();
           it2 != it1; ++it2)
        assert(comp(*it1, *it2) || comp(*it2, *it1));
    #endif
    std::sort(toOrder.begin(), toOrder.end(), comp);
  }

  // The nodes are collected, order them
  for (auto it = toOrder.begin(); it != toOrder.end(); ++it) {
    const Node *nd2 = *it;
    assert(nd1->getProcessID() != nd2->getProcessID() &&
           (isWrite(nd2) ||
            (isRead(nd2) && annotation.defines(nd2))) &&
           (!isRead(nd1) || !isRead(nd2)));

    // Prepare worklist
    worklist_ready.swap(worklist_done);
    assert(worklist_done.empty());

    if (worklist_ready.empty()) {
      // So far we work only with one po - the argument
      if (!areOrdered(nd1, nd2, po) &&
          (isRead(nd1) || isRead(nd2) ||
           (isObservable(nd1, po) && isObservable(nd2, po))
          )) {
        PartialOrder current = PartialOrder
          (std::unique_ptr<ThreadPairsVclocks>(new ThreadPairsVclocks(*(po.first))),
           std::unique_ptr<ThreadPairsVclocks>(new ThreadPairsVclocks(*(po.second))));
        PartialOrder otherorder = PartialOrder
          (std::unique_ptr<ThreadPairsVclocks>(new ThreadPairsVclocks(*(po.first))),
           std::unique_ptr<ThreadPairsVclocks>(new ThreadPairsVclocks(*(po.second))));
        // Handle current: w1 -> w2
        addEdge(nd1, nd2, current);
        worklist_done.push_back(std::move(current));
        assert(!current.first.get() &&
               !current.second.get());
        // Handle otherorder: w2 -> w1
        addEdge(nd2, nd1, otherorder);
        worklist_done.push_back(std::move(otherorder));
        assert(!otherorder.first.get() &&
               !otherorder.second.get());
      }
    } else {
      while (!worklist_ready.empty()) {
        PartialOrder current = std::move(worklist_ready.front());
        assert(!worklist_ready.front().first.get() &&
               !worklist_ready.front().second.get());
        worklist_ready.pop_front();

        if (!areOrdered(nd1, nd2, current) &&
            (isRead(nd1) || isRead(nd2) ||
             (isObservable(nd1, current) && isObservable(nd2, current))
            )) {
          // Unordered and either one of them read,
          // or both writes that are observable
          // Order these two
          PartialOrder otherorder = PartialOrder
            (std::unique_ptr<ThreadPairsVclocks>(new ThreadPairsVclocks(*(current.first))),
             std::unique_ptr<ThreadPairsVclocks>(new ThreadPairsVclocks(*(current.second))));
          // Handle current: w1 -> w2
          addEdge(nd1, nd2, current);
          worklist_done.push_back(std::move(current));
          assert(!current.first.get() &&
                 !current.second.get());
          // Handle otherorder: w2 -> w1
          addEdge(nd2, nd1, otherorder);
          worklist_done.push_back(std::move(otherorder));
          assert(!otherorder.first.get() &&
                 !otherorder.second.get());
        } else {
          // These two do not have to be ordered
          worklist_done.push_back(std::move(current));
          assert(!current.first.get() &&
                 !current.second.get());
        }
      }
    }
  }
}

std::vector<ZEvent> ZGraph::linearize(const PartialOrder& po,
                                              const ZAnnotation& annotation) const
{
  auto result = std::vector<ZEvent>();
  result.reserve(nodes_size());

  auto current = std::vector<unsigned>(processes.size(), 0);
  auto until = std::vector<unsigned>(processes.size(), 0);
  for (unsigned i=0; i<processes.size(); ++i) {
    const Node * last_nd_i = processes[i][ processes[i].size() - 1];
    if ((isRead(last_nd_i) && !annotation.defines(last_nd_i)) ||
        (isLock(last_nd_i) && !annotation.isLastLock(last_nd_i))) {
      // These nodes should be 'rediscovered' by the
      // TraceBuilder as something to annotate
      // Also locks can't be force-replayed if the
      // location is already held
      until[i] = processes[i].size() - 1;
    } else
      until[i] = processes[i].size();
  }
  ThreadPairsVclocks& pred = *(po.second);

  bool done = false;
  while (!done) {
    // Star-root process makes one step,
    // and before that step, all steps
    // of other processes that have
    // to HB because of 'po' are made
    auto requirements = std::map<unsigned, unsigned>();
    if (current[starRoot()] == until[starRoot()]) {
      // Star-root process finished, let all other processes finish
      for (unsigned i=0; i<processes.size(); ++i)
        if (i != starRoot() && current[i] < until[i])
          requirements.emplace(i, until[i]);
    } else {
      // Star-root process has not finished yet
      // Collect info of what needs to HB the next star-root step
      for (unsigned i=0; i<processes.size(); ++i)
        if (i != starRoot()) {
          int ipred = pred[starRoot()][i][ current[starRoot()] ] + 1;
          if (ipred > (int) current[i]) {
            requirements.emplace(i, ipred);
          }
        }
    }

    auto reqNodes = std::unordered_set<const Node *>();
    for (auto& rq : requirements) {
      assert(current[rq.first] < rq.second);
      assert(rq.second <= until[rq.first]);
      reqNodes.insert(processes[rq.first][ current[rq.first] ]);
    }
    while (!reqNodes.empty()) {
      // Find (an arbitrary) head of reqNodes
      auto it = reqNodes.begin();
      assert(it != reqNodes.end());
      const Node *headnd = *it;
      ++it;
      while (it != reqNodes.end()) {
        if (hasEdge(*it, headnd, po))
          headnd = *it;
        ++it;
      }
      // Perform the head, replace in reqNodes with its
      // thread successor (if that one is also required)
      result.push_back(headnd->getEvent()->copy(result.size(), true));
      unsigned tid = headnd->getProcessID();
      ++current[tid];
      assert(reqNodes.count(headnd));
      reqNodes.erase(headnd);
      if (current[tid] < requirements[tid])
        reqNodes.insert(processes[tid][ current[tid] ]);
    }

    if (current[starRoot()] < until[starRoot()]) {
      // Perform one step of the star-root process
      const Node * rootnd = processes[starRoot()][ current[starRoot()] ];
      result.push_back(rootnd->getEvent()->copy(result.size(), true));
      ++current[starRoot()];
    } else
      done = true;
  }

  return result;
}
*/
