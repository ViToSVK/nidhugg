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
#include <stack>

#include "VCHelpers.h"
#include "VCGraphVclock.h"
#include "SCC.h"

/* *************************** */
/* GRAPH EXTENSION             */
/* *************************** */

// Extends this graph so it corresponds to 'trace'
// Checks the header file for the method description
void VCGraphVclock::extendGraph(const std::vector<VCEvent>& trace) {
	
	assert(edges_worklist.empty());
	assert(!trace.empty());
	assert(trace.size() >= nodes_size());

	cpid_to_processid.reserve(trace.size() / 2);
	processes.reserve(trace.size() / 2);
	event_to_node.reserve(trace.size());
	read_vciid_to_node.reserve(trace.size());
	
	// Faster than cpid_to_processid
	std::unordered_map<int, unsigned> ipid_mapping;
  ipid_mapping.reserve(trace.size() / 2);
	
	std::vector<unsigned> cur_evidx(processes.size(), 0);
	cur_evidx.reserve(trace.size() / 2);
	
  std::vector<Node *> spawns;
	spawns.reserve(8);
  std::vector<Node *> joins;
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

	extension_from = INT_MAX;
  int trace_idx = -1;
  for (auto traceit = trace.begin(); traceit != trace.end(); ++traceit) {
		trace_idx++;
		const VCEvent *ev = &(*traceit);
		
		unsigned proc_idx = INT_MAX;

		auto ipidit = ipid_mapping.find(ev->iid.get_pid());
    if (ipidit != ipid_mapping.end()) {
      proc_idx = ipidit->second;
		} else {
      auto cpidit = cpid_to_processid.find(ev->cpid);
			if (cpidit != cpid_to_processid.end()) {
        proc_idx = cpidit->second;
				// add to ipid cache for faster lookup next time
				ipid_mapping.emplace(ev->iid.get_pid(), proc_idx);
			}
		}

		if (proc_idx == INT_MAX) {
      // A new process appeared
			proc_idx = processes.size();
			assert(cur_evidx.size() == proc_idx);
			cur_evidx.push_back(0);
			
      processes.emplace_back(std::vector<Node *>());
			processes[proc_idx].reserve(std::min((int) trace.size() / 2,
																					 (int) trace.size() - trace_idx));
			ipid_mapping.emplace(ev->iid.get_pid(), proc_idx);
			cpid_to_processid.emplace(ev->cpid, proc_idx);
		}
		
    unsigned ev_idx = cur_evidx[proc_idx];
    Node *nd = nullptr;
		
		if (ev_idx == processes[proc_idx].size()) {
			// New event
			if (extension_from == INT_MAX)
				extension_from = (unsigned) trace_idx;
			nd = new Node(proc_idx, ev_idx, ev);
			
			nodes.insert(nd);
			processes[proc_idx].push_back(nd);
			assert(processes[proc_idx].size() == ev_idx + 1);
			event_to_node.emplace(ev, nd);
			if (isRead(ev))
				read_vciid_to_node.emplace(VCIID(*ev), nd);
			
      // XXX: function pointer calls not handled
			if (isSpawn(ev))
				spawns.push_back(nd);
			if (isJoin(ev))
				joins.push_back(nd);

		} else {
			// Already known event
			assert(ev_idx < processes[proc_idx].size());
			nd = processes[proc_idx][ev_idx];
			assert(nd);
			assert(nodes.count(nd));
			assert(nd->getEvent.equalVCIID(&ev);
						 && "Original part of the basis is changed in 'trace'");
			// since VCIID's are equal, read_vciid_to_node is
			// already set up for 'nd' in case 'nd' is a read node

			auto etnit = event_to_node.find(nd->getEvent());
			assert(etnit != event_to_node.end());
			event_to_node.erase(etnit);
			event_to_node.emplace(ev, nd);

			nd->setEvent(ev);
		}

    ++cur_evidx[proc_idx];

		if (isRead(ev)) {
      if (!tw_candidate.count(ev->ml))
				tw_candidate.emplace(ev->ml, std::vector<std::vector<int>>());
		}
		if (isMutexInit(ev)) {
			assert(!mutex_inits.count(ev->ml));
			mutex_inits.emplace(ev->ml, nd);
		}
		if (isMutexDestroy(ev)) {
			assert(!mutex_destroys.count(ev->ml));
			mutex_destroys.emplace(ev->ml, nd);
		}
		if (isLock(ev) || isUnlock(ev)) {
			// for each ml, for each thread, remember first access
      if (!mutex_first.count(ev->ml)) {
				mutex_first.emplace(ev->ml,
														std::unordered_map<unsigned, const Node *>());
				mutex_first[ev->ml].reserve(8);
			}
      if (!mutex_first[ev->ml].count(nd->getProcessID()))
				mutex_first[ev->ml].emplace(nd->getProcessID(),
																		nd);
			// for each ml, for each thread, remember last access
      if (!mutex_last.count(ev->ml)) {
				mutex_last.emplace(ev->ml,
														std::unordered_map<unsigned, const Node *>());
				mutex_last[ev->ml].reserve(8);
			}
			// creates if doesn't exist, replaces if exists
      mutex_last[ev->ml][nd->getProcessID()] = nd;
		}

	}

	assert(processes.size() == cur_evidx.size());
	for (unsigned i = 0; i < cur_evidx.size(); ++i) {
    assert(curevidx[i] == processes[i].size()
					 && "Didn't go through entire original part of the basis");
	}
	
	// EDGES - extend for original processes

	succ_original.reserve(processes.size());
	pred_original.reserve(processes.size());	

  for (unsigned i=0; i<succ_original.size(); i++) {
		succ_original[i].reserve(processes.size());
		pred_original[i].reserve(processes.size());
    for (unsigned j=0; j<processes.size(); j++) {
			if (j < succ_original.size() && i != j) {
        // Vclocks from original to original
				succ_original[i][j].reserve(processes[i].size());
				pred_original[i][j].reserve(processes[i].size());
				// new succ slots should be filled with INT_MAX,
				// new pred slots should be filled with what the last original slot says
				int last_pred = pred_original[i][j]
					                           [pred_original[i][j].size() - 1];
				while (succ_original[i][j].size() != processes[i].size()) {
          succ_original[i][j].push_back(INT_MAX);
					pred_original[i][j].push_back(last_pred);
				}
				assert(succ_original[i][j].size() == processes[i].size() &&
							 pred_original[i][j].size() == processes[i].size());
			}
		  if (j >= succ_original.size()) {
        // Vclocks from original to new
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
	
	// EDGES - spawns
	for (const Node *spwn : spawns) {
		auto spwn_it = cpid_to_processid.find(spwn->getEvent()->childs_cpid);
		assert (spwn_it != cpid_to_processid.end());
		const Node *nthr = processes[spwn_it->second][0];
		
		assert(!hasEdge(nthr, spwn, &succ_original));
		if (!hasEdge(spwn, nthr, &succ_original))
      addEdge(spwn, nthr, &succ_original, &pred_original);
	}

	// EDGES - joins
	for (const Node *jn : joins) {
		auto jn_it = cpid_to_processid.find(jn->getEvent()->childs_cpid);
		assert (jn_it != cpid_to_processid.end());
    const Node *wthr = processes[jn_it->second]
			                          [processes[jn_it->second].size() - 1];

		assert(!hasEdge(jn, wthr, &succ_original));
		if (!hasEdge(wthr, jn, &succ_original))
	    addEdge(wthr, jn, &succ_original, &pred_original);
	}

	// EDGES - mutex inits
	for (auto& in : mutex_inits) {
    const SymAddrSize& loc_init = in.first;
		const Node *nd_init = in.second;

		for (auto& tid_nd_first : mutex_first[loc_init]) {
			const Node *nd_first = tid_nd_first.second;

			assert(!hasEdge(nd_first, nd_init, &succ_original));
			if (!hasEdge(nd_init, nd_first, &succ_original))
				addEdge(nd_init, nd_first, &succ_original, &pred_original);
		}
	}

	// EDGES - mutex destroys
	for (auto& de : mutex_destroys) {
    const SymAddrSize& loc_destroy = de.first;
		const Node *nd_destroy = de.second;

		for (auto& tid_nd_last : mutex_last[loc_destroy]) {
			const Node *nd_last = tid_nd_last.second;

			assert(!hasEdge(nd_destroy, nd_last, &succ_original));
			if (!hasEdge(nd_last, nd_destroy, &succ_original))
				addEdge(nd_last, nd_destroy, &succ_original, &pred_original);
		}
	}

	// TAIL WRITE CANDIDATES CACHE
	// [ml][tid][evid] returns idx of first event of thread-tid writing to ml
	// starting from AND INCLUDING evid and going back - (evid, evid-1, .., 0)
	// returns -1 if there is no such write
	
  for (auto& ml_cache : tw_candidate) {
    ml_cache.second.reserve(processes.size());
		for (unsigned i=0; i<processes.size(); ++i)
			ml_cache.second.push_back(std::vector<int>(processes[i].size(), -1));
	}

	for (unsigned tid = 0; tid < processes.size(); ++tid)
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
		}
	
}

/* *************************** */
/* EDGE ADDITION               */
/* *************************** */

// Adds an edge between two unordered nodes
// This method maintains:
// 1) thread-pair-wise transitivity
// 2) complete transitivity
void VCGraphVclock::addEdge(const Node *n1, const Node *n2,
														ThreadPairsVclocks *succ,
														ThreadPairsVclocks *pred) {
	assert(n1 && n2 && "Do not have such node");
	assert(!areOrdered(n1, n2, succ));
	unsigned ti = n1->getProcessID();
	unsigned tj = n2->getProcessID();
	unsigned ti_evx = n1->getEventID();
	unsigned tj_evx = n2->getEventID();
	
	addEdgeHelp(ti, ti_evx, tj, tj_evx, succ, pred);

  // Maintenance of the complete transitivity
  // Collect nodes from different threads with edges:
	// 1) to   ti[ti_evx]
	// 2) from tj[tj_evx]
	
	std::set< std::pair<unsigned, unsigned>> before_tiEvent, after_tjEvent;
	for (unsigned k = 0; k<processes.size(); ++k) {
    if (k != ti && k != tj) {
      int maxbefore = (*pred)[ti][k][ti_evx];
			if (maxbefore >= 0)
				before_tiEvent.emplace(k, maxbefore);

			int minafter = (*succ)[tj][k][tj_evx];
			assert(minafter >= 0);
			if ((unsigned) minafter < processes[k].size())
				after_tjEvent.emplace(k, minafter);
		}
	}

  // Add edges between each of 1) and tj[tj_evx]
	for (auto& bef : before_tiEvent)
		addEdgeHelp(bef.first, bef.second,
								tj, tj_evx,
								succ, pred);

	// Add edges between ti[ti_evx] and each of 2)

	for (auto& aft : after_tjEvent)
		addEdgeHelp(ti, ti_evx,
								aft.first, aft.second,
								succ, pred);
	
	// Add edges between each of 1) and each of 2)
	// (if they belong to different threads)
	for (auto& bef : before_tiEvent)
		for (auto& aft : after_tjEvent)
			if (bef.first != aft.first)
        addEdgeHelp(bef.first, bef.second,
										aft.first, aft.second,
										succ, pred);
	
}

// Helper method for addEdge,
// maintains thread-pair-wise transitivity
void VCGraphVclock::addEdgeHelp(unsigned ti, unsigned ti_evx,
																unsigned tj, unsigned tj_evx,
																ThreadPairsVclocks *succ,
																ThreadPairsVclocks *pred) {
	assert( (*succ)[ti][tj][ti_evx] > (int) tj_evx && // ! ti[ti_evx] HB tj[tj_evx]
					(*succ)[tj][ti][tj_evx] > (int) ti_evx && // ! tj[tj_evx] HB ti[ti_evx]
					"Tried to add an edge between ordered nodes");
	assert( (*pred)[tj][ti][tj_evx] < (int) ti_evx && // ! ti[ti_evx] HB tj[tj_evx]
					(*pred)[ti][tj][ti_evx] < (int) tj_evx && // ! tj[tj_evx] HB ti[ti_evx]
				 && "Inconsistent succ/pred vector clocks");

	(*succ)[ti][tj][ti_evx] = (int) tj_evx;
	(*pred)[tj][ti][tj_evx] = (int) ti_evx;
	
  // Maintenance of thread-pair-wise transitivity of succ[ti][tj]
	// Everything in ti before ti_evx also happens-before tj[tj_evx]
	
	for (int ti_before_evx = ti_evx - 1;
			 ti_before_evx >= 0;
			 --ti_before_evx) {
		if ((*succ)[ti][tj][ti_before_evx] <= (int) tj_evx)
      break; // since for smaller indices also <= tj_evx
		(*succ)[ti][tj][ti_before_evx] = (int) tj_evx;
	}

  // Maintenance of thread-pair-wise transitivity of pred[tj][ti]
	// Everything in tj after tj_evx also happens-after ti[ti_evx]
	
	for (unsigned tj_after_evx = tj_evx + 1;
			 tj_after_evx < processes[tj].size();
			 ++tj_after_evx) {
		if ((*pred)[tj][ti][tj_after_evx] >= (int) ti_evx)
			break; // since for bigger indices also >= ti_evx
		(*pred)[tj][ti][tj_after_evx] = (int) ti_evx;
	}
	
}

const Node * VCGraphVclock::getTWcandidate(const Node *nd, unsigned thr_id,
																					 const ThreadPairsVclocks *succ) const {
	int ev_id = (nd->getProcessID() == thr_id) ? nd->getEventID()
						: (*succ)[nd->getProcessID()][thr_id][nd->getEventID()];
	if (ev_id >= (int) processes[thr_id].size())
		ev_id = processes[thr_id].size() - 1;
	else
		--ev_id;

	while (ev_id >= 0) {
		const VCEvent *ev = processes[thr_id][ev_id]->getEvent();
		if (isWrite(ev) && ev->ml == nd->getEvent()->ml)
			return processes[thr_id][ev_id];
		--ev_id;
	}

	return nullptr;
}

std::pair<std::unordered_set<const Node *>, std::unordered_set<const Node *>>
	VCGraphVclock::tailWrites(const Node *nd, const ThreadPairsVclocks *succ, int good) const {
	assert(nd && "Do not have such node");
	assert(nd != initial_node && "Asking getTWcandidate for the initial node");
	assert(isRead(nd->getEvent()));

  auto result = std::pair<std::unordered_set<const Node *>, std::unordered_set<const Node *>>();  
	
	auto mayBeOrdered = std::unordered_set<const Node *>();

	for (unsigned thr_id = 0; thr_id < processes.size(); ++thr_id) {
    const Node *candidate = getTWcandidate(nd, thr_id, succ);
		if (!candidate) {
      assert(isWrite(candidate->getEvent()));
			mayBeOrdered.emplace(candidate);
		}
	}

  if (mayBeOrdered.empty()) {
    // the initial node is the unique tail write
		if (good == 0)
			result.first.emplace(initial_node);
		else
			result.second.emplace(initial_node);
		return result;
	}
	
	auto trueTails = std::unordered_set<const Node *>(mayBeOrdered);
	assert(trueTails.size() == mayBeOrdered.size());

	for (auto it_nd1 = mayBeOrdered.begin(); it_nd1 != mayBeOrdered.end(); ++it_nd1) {
    auto it_nd2 = it_nd1;
		++it_nd2;
		while (it_nd2 != mayBeOrdered.end()) {
			if (hasEdge(*it_nd1, *it_nd2, succ))
				trueTails.erase(*it_nd1);
			if (hasEdge(*it_nd2, *it_nd1, succ))
				trueTails.erase(*it_nd2);
      ++it_nd2;
		}
	}

	assert(!trueTails.empty());

	for (auto& tail : trueTails) {
    if (tail->getEvent()->value == good)
			result.first.emplace(tail);
		else
			result.second.emplace(tail);
	}
  
	
	return result;
}

bool VCGraphVclock::valueClosure(const VCAnnotation& annot,
																 ThreadPairsVclocks *succ,
																 ThreadPairsVclocks *pred) {
	std::vector<int> done; // done on all bigger indices
	done.reserve(processes.size());
	for (unsigned tid=0; tid<processes.size(); ++tid)
		done.push_back(processes[tid].size() - 1);
	
	while (true) {
		auto maxReadCandidates = std::unordered_set<const Node *>();
		
    for (unsigned tid=0; tid<processes.size(); ++tid) {
      while (done[tid] > -1) {
        const Node *nd = processes[tid][done[tid]];
				if (isRead(nd->getEvent())) {
          maxReadCandidates.emplace(nd);
					break;
				}
				--done[tid];
			}
		}

		// finished scanning for maximal not-done read candidates

		if (maxReadCandidates.empty()) {
      for (unsigned tid=0; tid<processes.size(); ++tid)
				assert(done[tid] == -1);
			return true;
		}

		// there are some candidates, so one will be found

		const Node *maxRead = nullptr;

		while (!maxRead) {
			assert(!maxReadCandidates.empty());
      auto it_nd1 = maxReadCandidates.begin();
			auto it_nd2 = maxReadCandidates.begin();
			++it_nd2;

			while (it_nd2 != maxReadCandidates.end()) {
        if (hasEdge(*it_nd1, *it_nd2, succ))
          break;
				++it_nd2;
			}
			
			if (it_nd2 == maxReadCandidates.end()) {
        maxRead = *it_nd1;
				maxReadCandidates.clear();
			} else
				maxReadCandidates.erase(*it_nd1);
		}

		// have a maximal not-done read event
		
    int good = annot.getValue(*(maxRead->getEvent()));
		
		auto tailW = tailWrites(maxRead, succ, good);

		while (tailW.first.empty()) {
			// only bad tail writes
			// order all after maxRead
			// or return false if impossible
      assert(!tailW.second.empty());
			for (const Node *badtw : tailW.second) {
        if (hasEdge(badtw, maxRead, succ))
					return false;
				addEdge(maxRead, badtw, succ, pred);
			}

			tailW = tailWrites(maxRead, succ, good);
		}

		assert(!tailW.first.empty());
		--done[maxRead->getProcessID()];
			
	}
}
