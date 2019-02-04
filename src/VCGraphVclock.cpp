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

#include <unordered_map>
#include <map>
#include <stack>

#include "VCHelpers.h"
#include "VCGraphVclock.h"

/* *************************** */
/* GRAPH EXTENSION             */
/* *************************** */

// Extends this graph so it corresponds to 'trace'
// Check the header file for the method description
void VCGraphVclock::extendGraph(const std::vector<VCEvent>& trace,
                                const VCAnnotation *annotationPtr)
{
  assert(worklist_ready.empty());
  assert(worklist_done.empty());
  assert(!trace.empty());
  assert(trace.size() >= nodes_size());

  cpid_to_processid.reserve(4);
  processes.reserve(4);
  event_to_node.reserve(trace.size());

  // Faster than cpid_to_processid
  std::unordered_map<int, unsigned> ipid_to_processid;
  ipid_to_processid.reserve(4);

  std::vector<unsigned> cur_evidx(processes.size(), 0);
  cur_evidx.reserve(4);

  std::unordered_map<unsigned, const Node *> has_unannotated_read_or_lock;
  has_unannotated_read_or_lock.reserve(4);

  std::unordered_set<unsigned> has_second_unannot;
  has_second_unannot.reserve(4);

  std::unordered_set<SymAddrSize> found_last_lock_for_location;
  found_last_lock_for_location.reserve(4);

  std::unordered_set<unsigned> forbidden_processes_ipid;
  forbidden_processes_ipid.reserve(4);

  std::unordered_set<CPid> processes_created_within_po_cpid;
  processes_created_within_po_cpid.reserve(4);
  processes_created_within_po_cpid.insert(trace[0].cpid);

  std::unordered_set<unsigned> processes_with_event_we_dont_add;
  processes_with_event_we_dont_add.reserve(4);

  scores_writeno.reserve(4);
  scores_conflict.reserve(4);
  scores_valueconflict.reserve(4);

  std::vector<const Node *> spawns;
  spawns.reserve(8);

  std::vector<const Node *> joins;
  joins.reserve(8);

  std::unordered_map<SymAddrSize, const Node *> mutex_inits;
  mutex_inits.reserve(8);

  std::unordered_map<SymAddrSize, const Node *> mutex_destroys;
  mutex_destroys.reserve(8);

  std::unordered_map<SymAddrSize,
                     std::unordered_map<unsigned, const Node *>> mutex_first;
  mutex_first.reserve(8);

  std::unordered_map<SymAddrSize,
                     std::unordered_map<unsigned, const Node *>> mutex_last;
  mutex_first.reserve(8);

  for (auto traceit = trace.begin(); traceit != trace.end(); ++traceit) {
    const VCEvent *ev = &(*traceit);

    if (forbidden_processes_ipid.count(ev->iid.get_pid())) {
        // This is a forbidden process because it
        // happens after a so-far unannotated node
        continue;
    }

    unsigned proc_idx = INT_MAX;

    // Check if this process is already known
    auto ipidit = ipid_to_processid.find(ev->iid.get_pid());
    if (ipidit != ipid_to_processid.end()) {
      proc_idx = ipidit->second;
    } else {
      auto cpidit = cpid_to_processid.find(ev->cpid);
      if (cpidit != cpid_to_processid.end()) {
        proc_idx = cpidit->second;
        // add to ipid cache for faster lookup next time
        ipid_to_processid.emplace(ev->iid.get_pid(), proc_idx);
      }
    }

    if (proc_idx == INT_MAX) {
      // A new process appeared
      if (!processes_created_within_po_cpid.count(ev->cpid)) {
        // This is a forbidden process because it
        // is created after a so-far unannotated node
        forbidden_processes_ipid.insert(ev->iid.get_pid());
        continue;
      }
      proc_idx = processes.size();
      assert(cur_evidx.size() == proc_idx);
      cur_evidx.push_back(0);

      processes.emplace_back(std::vector<Node *>());
      processes[proc_idx].reserve((int) trace.size() / 2);
      ipid_to_processid.emplace(ev->iid.get_pid(), proc_idx);
      cpid_to_processid.emplace(ev->cpid, proc_idx);
    }

    unsigned ev_idx = cur_evidx[proc_idx];

    // Check if we already haven't added a process-predecessor
    if (processes_with_event_we_dont_add.count(proc_idx)) {
      assert(ev_idx == processes[proc_idx].size());
      continue;
    }
    // Check if proc_idx already has an unannotated read
    if (has_unannotated_read_or_lock.count(proc_idx)) {
      // This event happens after something we need to
      // annotate first, so we don't add it into the graph
      assert(ev_idx == processes[proc_idx].size());
      if (has_second_unannot.count(proc_idx))
        continue;
      processes_with_event_we_dont_add.insert(proc_idx);
      if (isRead(ev) || isLock(ev)) {
        has_second_unannot.insert(proc_idx);
        continue;
      }
      // Compute scores for this thread
      if (isWrite(ev)) {
        if (scores_writeno.count(proc_idx))
          scores_writeno[proc_idx] = scores_writeno[proc_idx] + 1;
        else
          scores_writeno.emplace(proc_idx, 1);
        for (auto& prid_nd : has_unannotated_read_or_lock)
          if (isRead(prid_nd.second) && prid_nd.first != proc_idx) {
            if (sameMl(ev, prid_nd.second->getEvent())) {
              // Hides write that conflicts with nd
              if (!scores_conflict.count(proc_idx)) {
                scores_conflict.emplace(proc_idx, std::unordered_set<unsigned>());
                scores_valueconflict.emplace(proc_idx, std::unordered_set<unsigned>());
              }
              scores_conflict[proc_idx].insert(prid_nd.first);
              if (ev->value == prid_nd.second->getEvent()->value) {
                // Same value written than what nd observes now
                scores_valueconflict[proc_idx].insert(prid_nd.first);
              }
            }
          }
      }
      continue;
    }
    // Joins of processes with an event we didn't include also
    // can not be included and neither any of their successors
    if (isJoin(ev)) {
      auto childs_proc_it = cpid_to_processid.find(ev->childs_cpid);
      if (childs_proc_it == cpid_to_processid.end() ||
          processes_with_event_we_dont_add.count(childs_proc_it->second)) {
        assert(ev_idx == processes[proc_idx].size());
        processes_with_event_we_dont_add.insert(proc_idx);
        continue;
      }
    }

    Node *nd = nullptr;
    if (ev_idx == processes[proc_idx].size()) {
      // New event
      nd = new Node(proc_idx, ev_idx, ev);

      nodes.insert(nd);
      processes[proc_idx].push_back(nd);
      assert(processes[proc_idx].size() == ev_idx + 1);
      event_to_node.emplace(ev, nd);

      // XXX: function pointer calls not handled
      if (isSpawn(ev))
        spawns.push_back(nd);
      if (isJoin(ev))
        joins.push_back(nd);

    } else {
      // Already known event
      assert(ev_idx < processes[proc_idx].size());
      nd = processes[proc_idx][ev_idx];
      assert(nd && nodes.count(nd));
      #ifndef NDEBUG
      bool evType = (nd->getEvent()->kind == ev->kind);
      assert(evType && "Original part of the basis is changed in 'trace'");
      #endif
      auto etnit = event_to_node.find(nd->getEvent());
      assert(etnit != event_to_node.end());
      event_to_node.erase(etnit);
      event_to_node.emplace(ev, nd);

      nd->setEvent(ev);
    }

    assert(nd->getProcessID() == proc_idx &&
           nd->getEventID() == ev_idx &&
           nd->getEvent() == ev);
    assert(nd->getEvent()->event_order == nd->getEventID());
    nd->getEvent()->setPID(nd->getProcessID());

    ++cur_evidx[proc_idx];

    // Handle Spawn
    if (isSpawn(ev)) {
      processes_created_within_po_cpid.insert(ev->childs_cpid);
    }
    // Handle Write
    if (isWrite(ev)) {
      if (nd->getProcessID() == starRoot()) {
        // Root write
        auto itml = wRoot.find(ev->ml);
        if (itml == wRoot.end()) {
          wRoot.emplace_hint(itml, ev->ml, std::vector<const Node*>());
          wRoot[ev->ml].reserve(8);
        }
        wRoot[ev->ml].push_back(nd);
      } else {
        // Nonroot write
        leafThreadsWithRorW.insert(nd->getProcessID());
        auto itml = wNonrootUnord.find(ev->ml);
        if (itml == wNonrootUnord.end()) {
          wNonrootUnord.emplace_hint(itml, ev->ml, std::unordered_set<const Node*>());
          wNonrootUnord[ev->ml].reserve(8);
        }
        wNonrootUnord[ev->ml].insert(nd);
      }
    }

    // Handle Read
    if (isRead(ev)) {
      // Check if nd is not annotated yet
      if (!annotationPtr || !(annotationPtr->defines(nd))) {
        // This will be the last node for the corresponding thread
        has_unannotated_read_or_lock.emplace(nd->getProcessID(), nd);
      }
      auto ittw = tw_candidate.find(ev->ml);
      if (ittw == tw_candidate.end())
        tw_candidate.emplace_hint(ittw, ev->ml, std::vector<std::vector<int>>());
      auto itroot = wRoot.find(ev->ml);
      if (itroot == wRoot.end()) {
        wRoot.emplace_hint(itroot, ev->ml, std::vector<const Node*>());
        wRoot[ev->ml].reserve(8);
      }
      auto itnonr = wNonrootUnord.find(ev->ml);
      if (itnonr == wNonrootUnord.end()) {
        wNonrootUnord.emplace_hint(itnonr, ev->ml, std::unordered_set<const Node*>());
        wNonrootUnord[ev->ml].reserve(8);
      }
      if (nd->getProcessID() == starRoot()) {
        // Root read
        auto itml = readsRoot.find(ev->ml);
        if (itml == readsRoot.end()) {
          readsRoot.emplace_hint(itml, ev->ml, std::unordered_set<const Node*>());
          readsRoot[ev->ml].reserve(8);
        }
        readsRoot[ev->ml].insert(nd);
      } else {
        // Nonroot read
        leafThreadsWithRorW.insert(nd->getProcessID());
        auto itml = readsNonroot.find(ev->ml);
        if (itml == readsNonroot.end()) {
          readsNonroot.emplace_hint(itml, ev->ml, std::unordered_set<const Node*>());
          readsNonroot[ev->ml].reserve(8);
        }
        readsNonroot[ev->ml].insert(nd);
      }
    }

    // Handle Mutex events
    if (isMutexInit(ev)) {
      assert(!mutex_inits.count(ev->ml));
      mutex_inits.emplace(ev->ml, nd);
    }
    if (isMutexDestroy(ev)) {
      assert(!mutex_destroys.count(ev->ml));
      mutex_destroys.emplace(ev->ml, nd);
    }
    if (isLock(ev)) {
      // Check if nd is not annotated yet
      if (!annotationPtr || !(annotationPtr->locationHasSomeLock(nd))) {
        // No annotated lock has held this location yet
        // This will be the last node for the corresponding thread
        has_unannotated_read_or_lock.emplace(nd->getProcessID(), nd);
      } else {
        // Some lock has already happened on this location
        if (annotationPtr->isLastLock(nd)) {
          // This is the last annotated lock for this location
          assert(!found_last_lock_for_location.count(nd->getEvent()->ml));
          found_last_lock_for_location.insert(nd->getEvent()->ml);
        } else if (found_last_lock_for_location.count(nd->getEvent()->ml)) {
          // This is a lock after the last annotated lock for this location
          // This will be the last node for the corresponding thread
          has_unannotated_read_or_lock.emplace(nd->getProcessID(), nd);
        }
      }
    }
    if (isLock(ev) || isUnlock(ev)) {
      // For each ml, for each thread, remember first access
      if (!mutex_first.count(ev->ml)) {
        mutex_first.emplace(ev->ml, std::unordered_map<unsigned, const Node *>());
        mutex_first[ev->ml].reserve(8);
      }
      if (!mutex_first[ev->ml].count(nd->getProcessID()))
        mutex_first[ev->ml].emplace(nd->getProcessID(), nd);
      // For each ml, for each thread, remember last access
      if (!mutex_last.count(ev->ml)) {
        mutex_last.emplace(ev->ml, std::unordered_map<unsigned, const Node *>());
        mutex_last[ev->ml].reserve(8);
      }
      mutex_last[ev->ml][nd->getProcessID()] = nd;
    }
  } // end of loop for traversing trace and creating nodes

  processes.shrink_to_fit();
  for (unsigned i=0; i<processes.size(); ++i)
    processes[i].shrink_to_fit();

  #ifndef NDEBUG
  assert(processes.size() == cur_evidx.size());
  for (unsigned i = 0; i < cur_evidx.size(); ++i) {
    assert(cur_evidx[i] == processes[i].size()
           && "Didn't go through entire original part of the basis");
  }
  #endif

  // TAIL WRITE CANDIDATES CACHE
  // [ml][tid][evid] returns idx of first event of thread-tid writing to ml
  // starting from AND INCLUDING evid and going back - (evid, evid-1, .., 0)
  // returns -1 if there is no such write

  for (auto& ml_cache : tw_candidate) {
    ml_cache.second.reserve(processes.size());
    for (unsigned i=0; i<processes.size(); ++i)
      ml_cache.second.push_back(std::vector<int>(processes[i].size(), -1));
  }
  for (unsigned tid = 0; tid < processes.size(); ++tid) {
    for (unsigned evid = 0; evid < processes[tid].size(); ++evid) {
      // update [ml][tid][evid] for all ml
      const VCEvent *ev = processes[tid][evid]->getEvent();
      if (!isWrite(ev) || !tw_candidate.count(ev->ml)) {
        // no update anywhere
        if (evid > 0)
          for (auto& ml_cache : tw_candidate)
            ml_cache.second[tid][evid] = ml_cache.second[tid][evid-1];
      } else {
        // update on ev->ml
        bool foundalready = false;
        for (auto& ml_cache : tw_candidate)
          if (!foundalready && ml_cache.first == ev->ml) {
            foundalready = true;
            ml_cache.second[tid][evid] = evid;
          } else {
            if (evid > 0)
              ml_cache.second[tid][evid] = ml_cache.second[tid][evid-1];
          }
      }
    } // end of loop for evid
  } // end of loop for tid

  ThreadPairsVclocks& succ_original = *(original.first);
  ThreadPairsVclocks& pred_original = *(original.second);
  succ_original.reserve(processes.size());
  pred_original.reserve(processes.size());

  // CLOSURE SAFE UNTIL EVENT ID - init
  bool newProcesses = (succ_original.size() < processes.size());
  if (newProcesses) {
    closureSafeUntil = std::vector<int>(processes.size(), -1);
  } else {
    closureSafeUntil = std::vector<int>();
    closureSafeUntil.reserve(succ_original.size());
    for (unsigned i=0; i<succ_original.size(); i++) {
      closureSafeUntil.push_back( processes[i].size() - 2);
    }
  }

  // EDGES - extend for original processes
  for (unsigned i=0; i<succ_original.size(); i++) {
    succ_original[i].reserve(processes.size());
    pred_original[i].reserve(processes.size());
    for (unsigned j=0; j<processes.size(); j++) {
      if (j < succ_original.size() && i != j) {
        // Vclocks from original processes to original processes
        succ_original[i][j].reserve(processes[i].size());
        pred_original[i][j].reserve(processes[i].size());
        // new succ slots should be filled with INT_MAX,
        // new pred slots should be filled with what the last original slot says
        int last_pred = pred_original[i][j]
                                     [pred_original[i][j].size() - 1];
        if (succ_original[i][j].size() < processes[i].size()) {
          // New node at process i, all the nodes of process j
          // after last_pred become unsafe wrt closure
          if (closureSafeUntil[j] > last_pred)
            closureSafeUntil[j] = last_pred;
        }
        while (succ_original[i][j].size() < processes[i].size()) {
          succ_original[i][j].push_back(INT_MAX);
          pred_original[i][j].push_back(last_pred);
        }
        assert(succ_original[i][j].size() == processes[i].size() &&
               pred_original[i][j].size() == processes[i].size());
      }
      if (j >= succ_original.size()) {
        // Vclocks from original processes to new processes
        succ_original[i].push_back(std::vector<int>(processes[i].size(), INT_MAX));
        pred_original[i].push_back(std::vector<int>(processes[i].size(), -1));
      }
    }
  }

  // EDGES - create for new proccesses
  assert(succ_original.size() == pred_original.size());
  while (succ_original.size() < processes.size()) {
    unsigned i = succ_original.size();

    succ_original.push_back(std::vector<std::vector<int>>());
    succ_original[i].reserve(processes.size());
    pred_original.push_back(std::vector<std::vector<int>>());
    pred_original[i].reserve(processes.size());

    for (unsigned j=0; j<processes.size(); j++) {
      if (i != j) {
        succ_original[i].push_back(std::vector<int>(processes[i].size(), INT_MAX));
        pred_original[i].push_back(std::vector<int>(processes[i].size(), -1));
      } else {
        succ_original[i].push_back(std::vector<int>());
        pred_original[i].push_back(std::vector<int>());
      }
    }
  }

  succ_original.shrink_to_fit();
  pred_original.shrink_to_fit();
  for (unsigned i=0; i<processes.size(); ++i) {
    succ_original[i].shrink_to_fit();
    pred_original[i].shrink_to_fit();
  }

  // EDGES - spawns
  for (const Node *spwn : spawns) {
    auto spwn_it = cpid_to_processid.find(spwn->getEvent()->childs_cpid);
    assert (spwn_it != cpid_to_processid.end());
    const Node *nthr = processes[spwn_it->second][0];

    assert(!hasEdge(nthr, spwn, original));
    if (!hasEdge(spwn, nthr, original))
      addEdge(spwn, nthr, original);
  }

  // EDGES - joins
  for (const Node *jn : joins) {
    auto jn_it = cpid_to_processid.find(jn->getEvent()->childs_cpid);
    assert (jn_it != cpid_to_processid.end());
    const Node *wthr = processes[jn_it->second]
                                [processes[jn_it->second].size() - 1];

    assert(!hasEdge(jn, wthr, original));
    if (!hasEdge(wthr, jn, original))
      addEdge(wthr, jn, original);
  }

  // EDGES - mutex inits
  for (auto& in : mutex_inits) {
    const SymAddrSize& loc_init = in.first;
    const Node *nd_init = in.second;

    for (auto& tid_nd_first : mutex_first[loc_init]) {
      const Node *nd_first = tid_nd_first.second;

      assert(!hasEdge(nd_first, nd_init, original));
      if (!hasEdge(nd_init, nd_first, original))
        addEdge(nd_init, nd_first, original);
    }
  }

  // EDGES - mutex destroys
  for (auto& de : mutex_destroys) {
    const SymAddrSize& loc_destroy = de.first;
    const Node *nd_destroy = de.second;

    for (auto& tid_nd_last : mutex_last[loc_destroy]) {
      const Node *nd_last = tid_nd_last.second;

      assert(!hasEdge(nd_destroy, nd_last, original));
      if (!hasEdge(nd_last, nd_destroy, original))
        addEdge(nd_last, nd_destroy, original);
    }
  }
}

/* *************************** */
/* EDGE ADDITION               */
/* *************************** */

// Adds an edge between two unordered nodes
// This method maintains:
// 1) thread-pair-wise transitivity
// 2) complete transitivity
void VCGraphVclock::addEdge(const Node *n1, const Node *n2, const PartialOrder& po) const
{
  assert(n1 && n2 && "Do not have such node");
  assert(!areOrdered(n1, n2, po));
  unsigned ti = n1->getProcessID();
  unsigned tj = n2->getProcessID();
  unsigned ti_evx = n1->getEventID();
  unsigned tj_evx = n2->getEventID();

  addEdgeHelp(ti, ti_evx, tj, tj_evx, po);

  // Maintenance of the complete transitivity
  // Collect nodes from different threads with edges:
  // 1) to   ti[ti_evx]
  // 2) from tj[tj_evx]

  ThreadPairsVclocks& succ = *(po.first);
  ThreadPairsVclocks& pred = *(po.second);

  std::set< std::pair<unsigned, unsigned>> before_tiEvent, after_tjEvent;
  for (unsigned k = 0; k<processes.size(); ++k) {
    if (k != ti && k != tj) {
      int maxbefore = pred[ti][k][ti_evx];
      if (maxbefore >= 0)
        before_tiEvent.emplace(k, maxbefore);

      int minafter = succ[tj][k][tj_evx];
      assert(minafter >= 0);
      if ((unsigned) minafter < processes[k].size())
        after_tjEvent.emplace(k, minafter);
    }
  }

  // TryAdd edges between each of 1) (n3) and tj[tj_evx] (n2)
  for (auto& bef : before_tiEvent) {
    const Node *n3 = processes[bef.first][bef.second];
    assert(!hasEdge(n2, n3, po) && "Cycle");
    if (!hasEdge(n3, n2, po))
      addEdgeHelp(bef.first, bef.second,
                  tj, tj_evx, po);
  }

  // TryAdd edges between ti[ti_evx] (n1) and each of 2) (n3)
  for (auto& aft : after_tjEvent) {
    const Node *n3 = processes[aft.first][aft.second];
    assert(!hasEdge(n3, n1, po) && "Cycle");
    if (!hasEdge(n1, n3, po))
      addEdgeHelp(ti, ti_evx,
                  aft.first, aft.second, po);
  }

  // TryAdd edges between each of 1) (n3) and each of 2) (n4)
  // (if they belong to different threads)
  for (auto& bef : before_tiEvent)
    for (auto& aft : after_tjEvent)
      if (bef.first != aft.first) {
        const Node *n3 = processes[bef.first][bef.second];
        const Node *n4 = processes[aft.first][aft.second];
        assert(!hasEdge(n4, n3, po) && "Cycle");
        if (!hasEdge(n3, n4, po))
          addEdgeHelp(bef.first, bef.second,
                      aft.first, aft.second, po);
      }

}

// Helper method for addEdge,
// maintains thread-pair-wise transitivity
void VCGraphVclock::addEdgeHelp(unsigned ti, unsigned ti_evx,
                                unsigned tj, unsigned tj_evx,
                                const PartialOrder& po) const
{
  ThreadPairsVclocks& succ = *(po.first);
  ThreadPairsVclocks& pred = *(po.second);
  assert( succ[ti][tj][ti_evx] > (int) tj_evx && // ! ti[ti_evx] HB tj[tj_evx]
          succ[tj][ti][tj_evx] > (int) ti_evx && // ! tj[tj_evx] HB ti[ti_evx]
          "Tried to add an edge between ordered nodes");
  assert( pred[tj][ti][tj_evx] < (int) ti_evx && // ! ti[ti_evx] HB tj[tj_evx]
          pred[ti][tj][ti_evx] < (int) tj_evx && // ! tj[tj_evx] HB ti[ti_evx]
          "Inconsistent succ/pred vector clocks");

  // CLOSURE SAFE UNTIL EVENT ID - update before updating succ+pred
  // Must be done here because also edges added to maintain
  // transitivity can affect closure safety
  for (unsigned i=0; i<succ.size(); i++)
    if (i != tj) {
      // Events in thread i up until 'newbound' already happened
      // before [tj][tj_evx]
      int newbound = pred[tj][i][tj_evx];
      if (closureSafeUntil[i] > newbound)
        closureSafeUntil[i] = newbound;
    } else {
      if (closureSafeUntil[tj] > (int) tj_evx - 1)
        closureSafeUntil[tj] = (int) tj_evx - 1;
    }

  succ[ti][tj][ti_evx] = (int) tj_evx;
  pred[tj][ti][tj_evx] = (int) ti_evx;

  // Maintenance of thread-pair-wise transitivity of succ[ti][tj]
  // Everything in ti before ti_evx also happens-before tj[tj_evx]

  for (int ti_before_evx = ti_evx - 1;
       ti_before_evx >= 0;
       --ti_before_evx) {
    if (succ[ti][tj][ti_before_evx] <= (int) tj_evx)
      break; // since for smaller indices also <= tj_evx
    succ[ti][tj][ti_before_evx] = (int) tj_evx;
  }

  // Maintenance of thread-pair-wise transitivity of pred[tj][ti]
  // Everything in tj after tj_evx also happens-after ti[ti_evx]

  for (unsigned tj_after_evx = tj_evx + 1;
       tj_after_evx < processes[tj].size();
       ++tj_after_evx) {
    if (pred[tj][ti][tj_after_evx] >= (int) ti_evx)
      break; // since for bigger indices also >= ti_evx
    pred[tj][ti][tj_after_evx] = (int) ti_evx;
  }

}

/* *************************** */
/* MAIN ALGORITHM              */
/* *************************** */

const Node * VCGraphVclock::getTailWcandidate(const Node *nd, unsigned thr_id,
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
VCGraphVclock::getTailWrites(const Node *nd, const PartialOrder& po) const
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
VCGraphVclock::getHeadWcandidate(const Node *nd, unsigned thr_id, const PartialOrder& po) const
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
VCGraphVclock::getHeadWrites(const Node *nd, const PartialOrder& po) const
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

bool VCGraphVclock::isObservable(const Node *nd, const PartialOrder& po) const
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

bool VCGraphVclock::isObservableBy(const Node *writend, const Node *readnd,
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

void VCGraphVclock::orderEventMaz(const VCEvent *ev1, const VCAnnotation& annotation,
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

std::map<std::pair<int, VCAnnotation::Loc>, VCAnnotation::Ann>
VCGraphVclock::getMutationCandidates(const PartialOrder& po,
                                     const VCAnnotationNeg& negative, const Node *readnd) const
{
  assert(nodes.count(readnd));
  assert(isRead(readnd));

  int nd_tid = readnd->getProcessID();
  int nd_evid = readnd->getEventID();
  assert(nd_evid == (int) processes[nd_tid].size() - 1);

  auto mutateWrites = std::unordered_set<const Node *>();
  auto mayBeCovered = std::unordered_set<const Node *>();
  // From the value (evid) and below, everything is covered from nd by some other write
  auto covered = std::vector<int>(processes.size(), -1);

  ThreadPairsVclocks& succ = *(po.first);
  ThreadPairsVclocks& pred = *(po.second);

  // Handle thread nd_tid
  for (int evid = nd_evid - 1; evid >= 0; --evid) {
    const Node *wrnd = processes[nd_tid][evid];
    if (isWrite(wrnd) && sameMl(wrnd, readnd)) {
      // Add to mayBeCovered
      mayBeCovered.insert(wrnd);
      // Update cover in other threads
      for (unsigned i = 0; i < processes.size(); ++i)
        if (nd_tid != (int) i) {
          int newcov = pred[nd_tid][i][evid];
          if (newcov > covered[i])
            covered[i] = newcov;
        }
      // Nodes before wrnd in nd_tid are covered by wrnd
      break;
    }
  }

  // Handle all other threads
  for (unsigned tid = 0; tid < processes.size(); ++tid) {
    if (nd_tid != (int) tid) {
      // Handle nodes unordered with nd
      int su = succ[nd_tid][tid][nd_evid];
      if (su == INT_MAX)
        su = processes[tid].size();
      int pr = pred[nd_tid][tid][nd_evid];
      // (su-1, su-2, ..., pr+2, pr+1) unordered with nd
      for (int evid = su-1; evid > pr; --evid) {
        const Node *wrnd = processes[tid][evid];
        if (isWrite(wrnd) && sameMl(wrnd, readnd))
          mutateWrites.insert(wrnd);
        // wrnd covers nothing since it is unordered with nd
      }

      // Handle nodes that happen before nd
      for (int evid = pr; evid >= 0; --evid) {
        const Node *wrnd = processes[tid][evid];
        if (isWrite(wrnd) && sameMl(wrnd, readnd)) {
          // Add to mayBeCovered
          mayBeCovered.insert(wrnd);
          // Update cover in other threads
          for (unsigned i = 0; i < processes.size(); ++i)
            if (tid != i) {
              int newcov = pred[tid][i][evid];
              if (newcov > covered[i])
                covered[i] = newcov;
            }
          // Nodes before wrnd in tid are covered by wrnd
          break;
        }
      }
    }
  }

  // Only take those that are not covered from nd by some other
  for (auto& wrnd : mayBeCovered)
    if (covered[wrnd->getProcessID()] < (int) wrnd->getEventID())
      mutateWrites.emplace(wrnd);

  bool considerInitEvent = mayBeCovered.empty();

  // Disregard those forbidden by the negative annotation
  if (negative.forbidsInitialEvent(readnd)) {
    considerInitEvent = false;

    for (auto it = mutateWrites.begin(); it != mutateWrites.end(); ) {
      const Node * writend = *it;
      assert(isWrite(writend));
      if (negative.forbids(readnd, writend))
        it = mutateWrites.erase(it);
      else
        ++it;
    }
  }

  auto result = std::map<std::pair<int, VCAnnotation::Loc>,
                         VCAnnotation::Ann>();

  // Have mutateWrites, create possible annotations
  while (!mutateWrites.empty()) {
    auto it = mutateWrites.begin();
    const Node * writend = *it;

    int value = writend->getEvent()->value;
    VCAnnotation::Loc loc = (starRoot() != readnd->getProcessID())
      ? VCAnnotation::Loc::ANY :
      (writend->getProcessID() == readnd->getProcessID()
       ? VCAnnotation::Loc::LOCAL : VCAnnotation::Loc::REMOTE);

    if (loc == VCAnnotation::Loc::LOCAL) {
      // Done with this Ann
      auto goodRemote = std::unordered_set<VCIID>();
      auto goodLocal = VCIID(writend->getProcessID(), writend->getEventID());
      assert(!result.count(std::pair<int, VCAnnotation::Loc>(value, loc)));
      result.emplace(std::pair<int, VCAnnotation::Loc>(value, loc),
                     VCAnnotation::Ann(value, loc, std::move(goodRemote), true, goodLocal));
      mutateWrites.erase(it);
      continue;
    }

    assert(loc != VCAnnotation::Loc::LOCAL);

    auto goodRemote = std::unordered_set<VCIID>();
    auto goodLocal = VCIID(31337, 47);

    if (writend->getProcessID() == readnd->getProcessID()) {
      assert(loc == VCAnnotation::Loc::ANY);
      goodLocal = VCIID(writend->getProcessID(), writend->getEventID());
    }
    else
      goodRemote.emplace(writend->getProcessID(), writend->getEventID());

    if (considerInitEvent && value == 0 && loc == VCAnnotation::Loc::ANY) {
      assert(goodLocal.first == 31337);
      considerInitEvent = false;
      goodLocal = VCIID(INT_MAX, INT_MAX);
    }

    it = mutateWrites.erase(it);
    while (it != mutateWrites.end()) {
      const Node * anothernd = *it;
      if (value != anothernd->getEvent()->value ||
          (loc == VCAnnotation::Loc::REMOTE &&
           anothernd->getProcessID() == readnd->getProcessID())) {
        ++it;
      } else {
        // Good value and acceptable location
        if (anothernd->getProcessID() == readnd->getProcessID()) {
          assert(goodLocal.first == 31337);
          goodLocal = VCIID(anothernd->getProcessID(), anothernd->getEventID());
        }
        else
          goodRemote.emplace(anothernd->getProcessID(), anothernd->getEventID());
        it = mutateWrites.erase(it);
      }
    }

    assert(!result.count(std::pair<int, VCAnnotation::Loc>(value, loc)));
    result.emplace(std::pair<int, VCAnnotation::Loc>(value, loc),
                   VCAnnotation::Ann(value, loc,
                                     std::move(goodRemote),
                                     (goodLocal.first != 31337),
                                     goodLocal));
  }

  if (considerInitEvent) {
    VCAnnotation::Loc loc = (starRoot() != readnd->getProcessID())
      ? VCAnnotation::Loc::ANY : VCAnnotation::Loc::LOCAL;
    assert(!result.count(std::pair<int, VCAnnotation::Loc>(0, loc)));
    result.emplace(std::pair<int, VCAnnotation::Loc>(0, loc),
                   VCAnnotation::Ann(0, loc,
                                     std::unordered_set<VCIID>(),
                                     true,
                                     {INT_MAX, INT_MAX}));
  }

  return result;
}

std::vector<VCEvent> VCGraphVclock::linearize(const PartialOrder& po,
                                              const VCAnnotation& annotation) const
{
  auto result = std::vector<VCEvent>();
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
