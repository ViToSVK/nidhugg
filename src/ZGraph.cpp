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
  : tso(true), basis(*this), po(this->basis), cache()
{
  assert(&(basis.graph) == this);
  assert(&(po.basis) == &basis);
  assert(empty());
}


// Initial
ZGraph::ZGraph(const std::vector<ZEvent>& trace, int star_root_index, bool tso)
  : tso(tso),
    basis(*this, star_root_index),
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
  : tso(oth.tso),
    basis(std::move(oth.basis)),
    po(std::move(oth.po)),
    cache(std::move(oth.cache))
{
  assert(oth.empty());
}


// Partial order that will be moved as original
// Trace and annotation that will extend this copy of the graph
ZGraph::ZGraph
(const ZGraph& oth,
 ZPartialOrder&& po,
 const std::vector<ZEvent>& trace,
 const ZAnnotation& annotation)
  : tso(oth.tso),
    basis(oth.basis, *this),
    po(std::move(po), this->basis),
    cache()
{
  assert(&(this->basis.graph) == this);
  assert(&(this->po.basis) == &(this->basis));
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
    if (tso && ipidit != ipid_to_thraux.end()) {
      thr_idx = ipidit->second.first;
      assert(ev->auxID() == ipidit->second.second);
    } else {
      thr_idx = basis.getThreadID(ev).first;
      // add to ipid cache for faster lookup next time
      ipid_to_thraux.emplace(ev->iid.get_pid(),
                             std::pair<unsigned, int>(thr_idx,
                                                      ev->auxID()));
    }

    ev->_thread_id = thr_idx;
    if (!tso && isWriteM(ev)) {
      // Handle aux for pso
      int aux_idx = basis.psoGetAux(ev);
      ev->_aux_id = aux_idx;
    }

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

    if (!cur_evidx.count(thraux)) {
      assert(ev->eventID() == 0);
      cur_evidx.emplace(thraux, 0);
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
    //
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
    //
    if (isSpawn(ev)) {
      proc_seq_within_po.insert(ev->childs_cpid.get_proc_seq());
      assert(basis.number_of_threads() > 0);
      assert(ev->fence);
    }

    // Handle Buffer-Write
    //
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
    //
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
    }

    // Handle Read
    //
    if (isRead(ev)) {
      // Check if nd is not annotated yet
      if (!annotationPtr || !(annotationPtr->defines(ev))) {
        // This will be the last node for the corresponding thread
        has_unannotated_read_or_lock.insert(ev->threadID());
      } else {
        assert(annotationPtr->defines(ev));
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
    //
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
    for (int aux : basis.auxes(thr_added.first)) {
      assert(!basis(thr_added.first, aux).empty());
      const ZEvent *nthr = basis.getEvent(thr_added.first, aux, 0);

      assert(!po.hasEdge(nthr, spwn));
      if (!po.hasEdge(spwn, nthr))
        po.addEdge(spwn, nthr);
    }
  }

  // EDGES - joins
  for (const ZEvent *jn : joins) {
    auto thr_joined = basis.getThreadID(jn->childs_cpid.get_proc_seq());
    assert(!thr_joined.second);
    for (int aux : basis.auxes(thr_joined.first)) {
      assert(!basis(thr_joined.first, aux).empty());
      unsigned lastev_idx = basis(thr_joined.first, aux).size() - 1;
      const ZEvent *wthr = basis.getEvent(thr_joined.first, aux, lastev_idx);

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
}


/* *************************** */
/* MAIN ALGORITHM              */
/* *************************** */

int ZGraph::getTailWindex(const SymAddrSize& ml, unsigned thr, int evX) const
{
  assert(thr < basis.number_of_threads());
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
(const SymAddrSize& ml, unsigned thr, int evX) const
{
  int idx = getTailWindex(ml, thr, evX);
  if (idx == -1)
    return nullptr;

  assert(idx < (int) cache.wm.at(ml).at(thr).size());
  auto res = cache.wm.at(ml).at(thr)[idx];
  assert(isWriteM(res) && res->ml == ml &&
         res->threadID() == thr && res->auxID() != -1 &&
         evX >= 0 && res->eventID() <= (unsigned) evX);
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
      // Dumb
      /*
      int curr = su - 1;
      while (curr >= 0) {
        const ZEvent *rem = basis(tid, auxid)[curr];
        assert(isWriteM(rem) && !po.hasEdge(read, rem));
        if (po.hasEdge(rem, read)) {
          if (sameMl(read, rem)) {
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
          }
          curr--;
        } else {
          if (sameMl(read, rem))
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
    const ZEvent *localM = localB->write_other_ptr;
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
    // Delete remotes that happen *before* localM
    // (no chance for them, those after localM are safe)
    assert(isWriteM(localM) && sameMl(localB, localM) &&
           po.hasEdge(localB, localM));
    for (auto it = notCovered.begin(); it != notCovered.end(); ) {
      const ZEvent *rem = *it;
      assert(isWriteM(rem) && !po.areOrdered(read, rem));
      if (po.hasEdge(rem, localM))
        it = notCovered.erase(it);
      else
        ++it;
    }
    for (auto it = mayBeCovered.begin(); it != mayBeCovered.end(); ) {
      const ZEvent *rem = *it;
      assert(isWriteM(rem) && po.hasEdge(rem, read));
      if (po.hasEdge(rem, localM))
        it = mayBeCovered.erase(it);
      else
        ++it;
    }
    // Check if not covered by some remote
    bool localCovered = false;
    for (const auto& rem : mayBeCovered) {
      assert(isWriteM(rem) && po.hasEdge(rem, read));
      if (po.hasEdge(localM, rem)) {
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


std::vector<ZEvent> ZGraph::linearizeTSO
(const ZPartialOrder& partial, const ZAnnotation& annotation) const
{
  assert(tso);
  assert(basis.auxes(basis.root()).size() <= 2 && "TSO");

  std::vector<ZEvent> result;
  unsigned resultsize = 0;

  std::unordered_map<std::pair<unsigned, int>, unsigned>
    current, until;

  for (unsigned thr=0; thr<basis.number_of_threads(); ++thr)
    for (int aux : basis.auxes(thr)) {
      assert(aux <= 0);
      current.emplace(std::pair<unsigned, int>(thr, aux), 0);
      assert(!basis(thr, aux).empty());
      auto lastEv = basis(thr, aux)[basis(thr, aux).size() - 1];
      if ((isRead(lastEv) && !annotation.defines(lastEv)) ||
          (isLock(lastEv) && !annotation.isLastLock(lastEv))) {
        // These events should be 'rediscovered' by the
        // TraceBuilder as something to annotate
        // Also locks can't be force-replayed if the
        // location is already held
        until.emplace(std::pair<unsigned, int>(thr, aux),
                      basis(thr, aux).size() - 1);
        resultsize += basis(thr, aux).size() - 1;
      } else {
        until.emplace(std::pair<unsigned, int>(thr, aux),
                      basis(thr, aux).size());
        resultsize += basis(thr, aux).size();
      }
    }
  result.reserve(resultsize);

  assert(basis.hasThreadAux(basis.root(), -1));
  bool rootBufferExists =
    basis.hasThreadAux(basis.root(), 0);

  bool done = false;
  while (!done) {
    std::list<const ZEvent *> leaves, laux;

    for (unsigned thr=0; thr<basis.number_of_threads(); ++thr) {
      if (!basis.isRoot(thr)) {
        for (int aux : basis.auxes(thr)) {
          std::pair<unsigned, int> thau(thr, aux);
          if (current[thau] < until[thau]) {
            auto ev = basis(thau)[current[thau]];
            if (aux == -1)
              leaves.push_back(ev);
            else {
              assert(aux == 0);
              laux.push_front(ev);
            }
          }
        }
      }
    }
    for (const auto& la : laux)
      leaves.push_front(la);

    // Some leaf thread makes one step,
    // and before that step, all steps
    // of other thraux that have
    // to HB because of 'partial' are made

    const ZEvent *noreq = nullptr;
    const ZEvent *minreq = nullptr;
    int req = INT_MAX;

    for (const auto& lev : leaves) {
      // Is lev minimal from leaves?
      bool minimal = true;
      for (const auto& loth : leaves) {
        if (lev != loth && partial.hasEdge(loth, lev)) {
          minimal = false;
          break;
        }
      }

      if (!minimal)
        continue;

      assert(minimal);
      // Collect info what root buffer needs to HB this event
      if (!rootBufferExists) {
        // No root buffer events, we have a no-requirement event
        noreq = lev;
        break;
      } else {
        std::pair<unsigned, int> rootthau(basis.root(), 0);
        int id = partial.pred(lev, basis.root(), 0).second;
        assert(id + 1 <= (int) until[rootthau]);
        if (id < (int) current[rootthau]) {
          // We have a no-requirement event
          noreq = lev;
          break;
        } else if (id < req) {
          // We have a new minimal-requirement event
          minreq = lev;
          req = id;
        }
      }
    }

    if (noreq || !rootBufferExists) {
      // There are be only rootMain requirements
      std::pair<unsigned, int> rootthau(basis.root(), -1);
      int id = !noreq ? until[rootthau] - 1
        : partial.pred(noreq, basis.root(), -1).second;
      assert(id + 1 <= (int) until[rootthau]);
      while (id >= (int) current[rootthau]) {
        const ZEvent *rMain = basis(rootthau)[current[rootthau]];
        result.push_back(rMain->copy(result.size(), true));
        ++current[rootthau];
      }
      if (noreq) {
        result.push_back(noreq->copy(result.size(), true));
        std::pair<unsigned, int> leafthau(noreq->threadID(), noreq->auxID());
        assert(basis.hasThreadAux(leafthau) && current.count(leafthau));
        ++current[leafthau];
      } else {
        // Leaves are finished and there is no rootBuffer
        assert(!noreq && !minreq && !rootBufferExists && leaves.empty());
        assert(result.size() == resultsize);
        done = true;
      }
    } else {
      // There are rootMain and rootBuffer requirements
      assert(!noreq && rootBufferExists && (minreq || leaves.empty()));
      std::pair<unsigned, int> rMainThau(basis.root(), -1);
      std::pair<unsigned, int> rBufferThau(basis.root(), 0);
      int rMainId = !minreq ? until[rMainThau] - 1
        : partial.pred(minreq, basis.root(), -1).second;
      assert(rMainId + 1 <= (int) until[rMainThau]);
      int rBuffferId = !minreq ? until[rBufferThau] - 1
        : partial.pred(minreq, basis.root(), 0).second;
      assert(rBuffferId + 1 <= (int) until[rBufferThau]);
      while (rMainId >= (int) current[rMainThau] ||
             rBuffferId >= (int) current[rBufferThau]) {
        const ZEvent * rMain = rMainId >= (int) current[rMainThau]
          ? basis(rMainThau)[current[rMainThau]] : nullptr;
        const ZEvent * rBuffer = rBuffferId >= (int) current[rBufferThau]
          ? basis(rBufferThau)[current[rBufferThau]] : nullptr;
        assert(rMain || rBuffer);
        if (rBuffer && (!rMain || !partial.hasEdge(rMain, rBuffer))) {
          result.push_back(rBuffer->copy(result.size(), true));
          ++current[rBufferThau];
        } else {
          result.push_back(rMain->copy(result.size(), true));
          ++current[rMainThau];
        }
      }
      if (minreq) {
        result.push_back(minreq->copy(result.size(), true));
        std::pair<unsigned, int> leafthau(minreq->threadID(), minreq->auxID());
        assert(basis.hasThreadAux(leafthau) && current.count(leafthau));
        ++current[leafthau];
      } else {
        // Finished
        assert(result.size() == resultsize);
        done = true;
      }
    }
  }

  //dumpTrace(result);
  return result;
}


std::vector<ZEvent> ZGraph::linearizePSO
(const ZPartialOrder& partial, const ZAnnotation& annotation) const
{
  assert(!tso);
  ZPartialOrder lin = copyPO(partial);
  unsigned root = basis.root();

  // (1) For all conflicting unordered pairs (i) leaf-event (ii) root memory-write,
  //     order the leaf event before the root memory-write
  for (unsigned thr = 0; thr < basis.number_of_threads(); ++thr) {
    if (thr != root) {
      for (int aux : basis.auxes(thr)) {
        if (aux != -1) {
          // Leaf buffer
          assert(!basis(thr, aux).empty() && isWriteM(basis(thr, aux)[0]));
          const auto& ml = basis(thr, aux)[0]->ml;
          int raux = basis.auxForMl(ml, root);
          if (raux != -1) {
            // All events of thr,aux go before unordered events from root,raux
            for (int ev_id = basis(thr, aux).size() - 1; ev_id >= 0; --ev_id) {
              const ZEvent *leaf = basis(thr, aux)[ev_id];
              int pred = lin.pred(leaf, root, raux).second;
              if (pred + 1 < (int) basis(root, raux).size()) {
                const ZEvent *newsucc = basis(root, raux)[pred + 1];
                assert(!lin.hasEdge(newsucc, leaf));
                if (!lin.hasEdge(leaf, newsucc))
                  lin.addEdge(leaf, newsucc);
              }
            }
          }
        } else {
          // Leaf main thread
          for (int ev_id = basis(thr, aux).size() - 1; ev_id >= 0; --ev_id) {
            const ZEvent *leaf = basis(thr, aux)[ev_id];
            if (!isRead(leaf) || !annotation.defines(leaf)) // No need for buffer-writes
              continue;
            int raux = basis.auxForMl(leaf->ml, root);
            if (raux == -1)
              continue;
            int pred = lin.pred(leaf, root, raux).second;
            if (pred + 1 < (int) basis(root, raux).size()) {
              const ZEvent *newsucc = basis(root, raux)[pred + 1];
              assert(!lin.hasEdge(newsucc, leaf));
              if (!lin.hasEdge(leaf, newsucc))
                lin.addEdge(leaf, newsucc);
            }
          }
        }
      }
    }
  }

  // (2) Linearize arbitrarily
  std::vector<ZEvent> result;
  unsigned resultsize = 0;

  std::unordered_map<std::pair<unsigned, int>, unsigned>
    until;
  std::set<const ZEvent *, ZEventPtrComp> next;

  for (unsigned thr=0; thr<basis.number_of_threads(); ++thr)
    for (int aux : basis.auxes(thr)) {
      assert(!basis(thr, aux).empty());
      auto lastEv = basis(thr, aux)[basis(thr, aux).size() - 1];
      if ((isRead(lastEv) && !annotation.defines(lastEv)) ||
          (isLock(lastEv) && !annotation.isLastLock(lastEv))) {
        // These events should be 'rediscovered' by the
        // TraceBuilder as something to annotate
        // Also locks can't be force-replayed if the
        // location is already held
        until.emplace(std::pair<unsigned, int>(thr, aux),
                      basis(thr, aux).size() - 1);
        resultsize += basis(thr, aux).size() - 1;
        if (((int) basis(thr, aux).size()) - 1 > 0)
          next.insert(basis(thr, aux)[0]);
      } else {
        until.emplace(std::pair<unsigned, int>(thr, aux),
                      basis(thr, aux).size());
        resultsize += basis(thr, aux).size();
        next.insert(basis(thr, aux)[0]);
      }
    }
  result.reserve(resultsize);

  while (!next.empty()) {
    // Find an arbitrary head event
    const ZEvent *head = nullptr;
    for (const auto& ev : next) {
      if (!head || lin.hasEdge(ev, head))
        head = ev;
    }
    #ifndef NDEBUG
    for (const auto& ev : next)
      assert(ev == head || !lin.hasEdge(ev, head));
    #endif

    assert(next.count(head));
    next.erase(head);
    result.push_back(head->copy(result.size(), true));
    std::pair<unsigned, int> thau(head->threadID(), head->auxID());
    assert(until.count(thau));
    if (head->eventID() + 1 < until[thau]) {
      assert(head->eventID() + 1 < basis(thau).size());
      next.insert(basis(thau)[head->eventID() + 1]);
    }
  }

  assert(result.size() == resultsize);
  //dumpTrace(result);
  return result;
}


bool ZGraph::isObservable(const ZEvent *ev, const ZPartialOrder& partial) const
{
  assert(isWriteM(ev));
  assert(isWriteB(ev->write_other_ptr));

  // In local thread
  for (unsigned evIdx = ev->write_other_ptr->eventID() + 1;
       evIdx < basis(ev->threadID(), -1).size(); ++evIdx) {
    const ZEvent *read = basis(ev->threadID(), -1)[evIdx];
    if (isWriteB(read) && sameMl(ev, read)) {
      // no read between ev and the subsequent conflicting write
      break;
    }
    if (!isRead(read) || !sameMl(ev, read))
      continue;
    assert(cache.readWB.at(read) == ev->write_other_ptr);
    // ev observable by read if unordered or happening after read
    if (!partial.hasEdge(ev, read))
      return true;
    // might be covered from read by write in other thread
    assert(isRead(read) && sameMl(ev, read) && partial.hasEdge(ev, read));
    bool evCovered = isCoveredFromByOtherThread(ev, read, partial);
    if (!evCovered) {
      // ev observable by read
      return true;
    } else {
      // ev covered from read and hence covered
      // also from all potential local reads after
      // this one, stop searching in local thread
      break;
    }
  }

  // In remote threads
  for (unsigned thIdx = 0; thIdx < basis.number_of_threads(); ++thIdx) {
    if (thIdx == ev->threadID())
      continue;

    // Search for conflicting reads in thIdx
    for (int evIdx = partial.pred(ev, thIdx, -1).second + 1;
         evIdx < (int) basis(thIdx, -1).size(); ++evIdx) {
      const ZEvent *read = basis(thIdx, -1)[evIdx];
      assert(!partial.hasEdge(read, ev));
      if (!isRead(read) || !sameMl(ev, read))
        continue;

      if (!partial.hasEdge(ev, read)) {
        // Unordered, ev can be covered from read
        // only by a read-local write
        assert(!partial.areOrdered(ev, read));
        const ZEvent *writeRlocal = cache.readWB.at(read);
        if (!writeRlocal) {
          // ev observable by read
          return true;
        }
        assert(isWriteB(writeRlocal));
        writeRlocal = writeRlocal->write_other_ptr;
        assert(isWriteM(writeRlocal));
        if (!partial.hasEdge(ev, writeRlocal) &&
            !partial.hasEdge(read, writeRlocal)) {
          // ev observable by read
          return true;
        }
        // Not observable
        if (partial.hasEdge(read, writeRlocal)) {
          // ev not observable by read, but could be
          // by some other read in this same thread,
          // continue searching
          continue;
        } else {
          assert(partial.hasEdge(ev, writeRlocal));
          // since ev -> writeRlocal, it will be covered
          // by writeRlocal also to all later reads
          // in this thread
          break;
        }

      } else {
        // ev -> read, last chance in this thread
        // (1) Not covered by other threads
        if (isCoveredFromByOtherThread(ev, read, partial)) {
          // since ev -> read, it has no chance with reads after read in this thread
          break;
        }
        // (2) Not covered by ev-thread
        int auxIdx = basis.auxForMl(read->ml, ev->threadID());
        assert(auxIdx != -1 && auxIdx == ev->auxID());
        int pr = partial.pred(read, ev->threadID(), auxIdx).second;
        const ZEvent *othW = getTailW(read->ml, ev->threadID(), pr);
        assert(othW && isWriteM(othW));
        if (othW != ev) {
          // ev -> othW -> read
          assert(sameMl(othW, ev) && othW->eventID() > ev->eventID() &&
                 partial.hasEdge(othW, read));
          // since ev -> read, it has no chance with reads after read in this thread
          break;
        }
        // (3) Not covered by read-thread
        const ZEvent *writeRlocal = cache.readWB.at(read);
        if (!writeRlocal) {
          // ev observable by read
          return true;
        }
        assert(isWriteB(writeRlocal));
        writeRlocal = writeRlocal->write_other_ptr;
        assert(isWriteM(writeRlocal));
        if (!partial.hasEdge(ev, writeRlocal) &&
            !partial.hasEdge(read, writeRlocal)) {
          // ev observable by read
          return true;
        }
        assert(partial.hasEdge(ev, writeRlocal) ||
               partial.hasEdge(read, writeRlocal));
        // since ev -> read, it has no chance with reads after read in this thread
        break;
      }

    } // evIdx loop
  } // thIdx loop

  return false;
}


// Is 'writeEv' covered from 'readEv' in 'partial'
// by some other conflicting write othW? (writeEv -> othW -> readEv)
// This function checks for othW only in threads different from
// thread of writeEv and thread of readEv
bool ZGraph::isCoveredFromByOtherThread
(const ZEvent *writeEv, const ZEvent *readEv,
 const ZPartialOrder& partial) const
{
  assert(isWriteM(writeEv) && isRead(readEv) &&
         sameMl(writeEv, readEv) && partial.hasEdge(writeEv, readEv));
  for (unsigned thIdx = 0; thIdx < basis.number_of_threads(); ++thIdx) {
    // Only checks different threads!!!
    if (thIdx == writeEv->threadID() || thIdx == readEv->threadID())
      continue;
    int auxIdx = basis.auxForMl(readEv->ml, thIdx);
    if (auxIdx == -1)
      continue;
    int pr = partial.pred(readEv, thIdx, auxIdx).second;
    const ZEvent *othW = getTailW(readEv->ml, thIdx, pr);
    if (!othW)
      continue;
    assert(isWriteM(othW));
    assert(sameMl(writeEv, othW));
    assert(partial.hasEdge(othW, readEv));
    if (partial.hasEdge(writeEv, othW)) {
      // writeEv -> othW -> readEv
      // Hence writeEv covered
      return true;
    }
  }
  return false;
}
