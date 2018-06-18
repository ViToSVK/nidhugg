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
void VCGraphVclock::extendGraph(const std::vector<VCEvent>& trace)
{
	assert(worklist_ready.empty());
	assert(worklist_done.empty());
	assert(!trace.empty());
	assert(trace.size() >= nodes_size());

	cpid_to_processid.reserve(trace.size() / 2);
	processes.reserve(trace.size() / 2);
	event_to_node.reserve(trace.size());
	lock_vciid_to_node.reserve(trace.size());
	
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
			if (isLock(ev))
				lock_vciid_to_node.emplace(VCIID(*ev), nd);
			
      // XXX: function pointer calls not handled
			if (isSpawn(ev))
				spawns.push_back(nd);
			if (isJoin(ev))
				joins.push_back(nd);

		} else {
			// Already known event
			assert(extension_from == INT_MAX);
			assert(ev_idx < processes[proc_idx].size());
			nd = processes[proc_idx][ev_idx];
			assert(nd);
			assert(nodes.count(nd));
			assert(nd->getEvent()->equalVCIID(*ev)
						 && "Original part of the basis is changed in 'trace'");
			// since VCIID's are equal, lock_vciid_to_node is
			// already set up for 'nd' in case 'nd' is a lock node

			auto etnit = event_to_node.find(nd->getEvent());
			assert(etnit != event_to_node.end());
			event_to_node.erase(etnit);
			event_to_node.emplace(ev, nd);

			nd->setEvent(ev);
		}

    ++cur_evidx[proc_idx];

		if (isWrite(ev) && nd->getProcessID() != starRoot()) {
      auto itml = nonstar_writes.find(ev->ml);
			if (itml == nonstar_writes.end())
				nonstar_writes.emplace(ev->ml, std::unordered_set<const Node*>());
			nonstar_writes[ev->ml].insert(nd);
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
    assert(cur_evidx[i] == processes[i].size()
					 && "Didn't go through entire original part of the basis");
	}
	
	// EDGES - extend for original processes

  ThreadPairsVclocks& succ_original = *(original.first);
	ThreadPairsVclocks& pred_original = *(original.second);
	
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
void VCGraphVclock::addEdge(const Node *n1, const Node *n2, const PartialOrder& po)
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
				if (!hasEdge(n4, n3, po))
					addEdgeHelp(bef.first, bef.second,
											aft.first, aft.second, po);
			}

}

// Helper method for addEdge,
// maintains thread-pair-wise transitivity
void VCGraphVclock::addEdgeHelp(unsigned ti, unsigned ti_evx,
																unsigned tj, unsigned tj_evx,
																const PartialOrder& po)
{
  ThreadPairsVclocks& succ = *(po.first);
	ThreadPairsVclocks& pred = *(po.second);
	assert( succ[ti][tj][ti_evx] > (int) tj_evx && // ! ti[ti_evx] HB tj[tj_evx]
					succ[tj][ti][tj_evx] > (int) ti_evx && // ! tj[tj_evx] HB ti[ti_evx]
					"Tried to add an edge between ordered nodes");
	assert( pred[tj][ti][tj_evx] < (int) ti_evx && // ! ti[ti_evx] HB tj[tj_evx]
					pred[ti][tj][ti_evx] < (int) tj_evx && // ! tj[tj_evx] HB ti[ti_evx]
				  "Inconsistent succ/pred vector clocks");

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

bool VCGraphVclock::valueClosure(const VCAnnotation& annot,
																 const PartialOrder& po)
{
  return true;
	// TODO
}

void VCGraphVclock::orderEventMaz(const VCEvent *ev1)
{
	assert(isWrite(ev1) || isRead(ev1));
	auto it = event_to_node.find(ev1);
	assert(it != event_to_node.end());
	const Node *nd1 = it->second;
	assert(nd1->getEvent() == ev1);
	assert(nd1->getProcessID() != starRoot());

	auto itnsw = nonstar_writes.find(ev1->ml);
	if (itnsw == nonstar_writes.end())
		return;

	std::unordered_set<const Node *>& toOrder = itnsw->second;
	
	for (auto it = toOrder.begin(); it != toOrder.end(); ++it) {
    const Node *nd2 = *it;
		if (nd1 == nd2)
			continue;

		worklist_ready.swap(worklist_done);
		assert(worklist_done.empty());
		assert(!worklist_ready.empty());

		while (!worklist_ready.empty()) {
			PartialOrder current = std::move(worklist_ready.front());
			assert(!worklist_ready.front().first.get());
			assert(!worklist_ready.front().second.get());
			worklist_ready.pop_front();

      if (areOrdered(nd1, nd2, current)) {
        worklist_done.push_back(std::move(current));
				assert(!current.first.get());
				assert(!current.second.get());
			} else {
				PartialOrder otherorder = PartialOrder
					(std::unique_ptr<ThreadPairsVclocks>(new ThreadPairsVclocks(*(current.first))),
					 std::unique_ptr<ThreadPairsVclocks>(new ThreadPairsVclocks(*(current.second))));
				// handle current: w1 -> w2
				addEdge(nd1, nd2, current);
        worklist_done.push_back(std::move(current));
				assert(!current.first.get());
				assert(!current.second.get());
				// handle otherorder: w2 -> w1
				addEdge(nd2, nd1, otherorder);
        worklist_done.push_back(std::move(otherorder));
				assert(!otherorder.first.get());
				assert(!otherorder.second.get());
			}
		}
	}
}

std::unordered_set<int> VCGraphVclock::getMutateValues(const PartialOrder& po, const Node *nd) const
{	
	assert(nodes.count(nd));
	assert(isRead(nd->getEvent()));

	int nd_tid = nd->getProcessID();
	int nd_evid = nd->getEventID();
	assert(nd_evid == (int) processes[nd_tid].size() - 1);
	
	auto mutateWrites = std::unordered_set<const Node *>();
	auto mayBeCovered = std::unordered_set<const Node *>();
	// from here and below everything is covered from nd by some other write
	auto covered = std::vector<int>(processes.size(), -1);

  ThreadPairsVclocks& succ = *(po.first);
	ThreadPairsVclocks& pred = *(po.second);
	
	// handle thread nd_tid
	for (int evid = nd_evid - 1; evid >= 0; --evid) {
    const Node *wrnd = processes[nd_tid][evid];
		if (isWrite(wrnd->getEvent()) && wrnd->getEvent()->ml == nd->getEvent()->ml) {
			// add to mayBeCovered
			mayBeCovered.insert(wrnd);
      // update cover in other threads
			for (unsigned i = 0; i < processes.size(); ++i)
				if (nd_tid != (int) i) {
					int newcov = pred[nd_tid][i][evid];
					if (newcov > covered[i])
						covered[i] = newcov;
			  }
			// nodes before wrnd in nd_tid are covered by wrnd
			break;
		}
	}

	// handle all other threads
	for (unsigned tid = 0; tid < processes.size(); ++tid) {
		if (nd_tid != (int) tid) {
      // handle nodes unordered with nd
			int su = succ[nd_tid][tid][nd_evid];
			if (su == INT_MAX)
				su = processes[tid].size();
			int pr = pred[nd_tid][tid][nd_evid];
			// (su-1, su-2, ..., pr+2, pr+1) unordered with nd
			for (int evid = su-1; evid > pr; --evid) {
				const Node *wrnd = processes[tid][evid];
				if (isWrite(wrnd->getEvent()) && wrnd->getEvent()->ml == nd->getEvent()->ml)
					mutateWrites.insert(wrnd);
				// wrnd covers nothing since it is unordered with nd
			}

			// handle nodes that happen before nd
			for (int evid = pr; evid >= 0; --evid) {
				const Node *wrnd = processes[tid][evid];
				if (isWrite(wrnd->getEvent()) && wrnd->getEvent()->ml == nd->getEvent()->ml) {
					// add to mayBeCovered
					mayBeCovered.insert(wrnd);
					// update cover in other threads
					for (unsigned i = 0; i < processes.size(); ++i)
						if (tid != i) {
							int newcov = pred[tid][i][evid];
							if (newcov > covered[i])
								covered[i] = newcov;
						}
					// nodes before wrnd in tid are covered by wrnd
					break;
				}
			}
		}
	}

  // only take those that are not
	// covered from nd by some other
	for (auto& wrnd : mayBeCovered)
		if (covered[wrnd->getProcessID()] < (int) wrnd->getEventID())
			mutateWrites.emplace(wrnd);

	// have mutateWrites, return set of their values
  auto result = std::unordered_set<int>();
	for (auto& wrnd : mutateWrites)
		result.insert(wrnd->getEvent()->value);

	// if we have nothing to cover the init event
	if (mayBeCovered.empty())
		result.insert(0);
	
	return result;
}

std::vector<VCEvent> VCGraphVclock::linearize(const PartialOrder& po) const
{
	auto result = std::vector<VCEvent>();
	result.reserve(nodes_size());
	
	auto current = std::vector<unsigned>(processes.size(), 0);
	ThreadPairsVclocks& pred = *(po.second);

  bool done = false;	
	while (!done) {
		// star-root process makes one step,
		// and before that step, all steps
		// of other processes that have
		// to HB because of 'po' are made
		auto requirements = std::map<unsigned, unsigned>();
    if (current[starRoot()] == processes[starRoot()].size()) {
      // star-root process finished, let all other processes finish
			for (unsigned i=0; i<processes.size(); ++i)
				if (i != starRoot() && current[i] < processes[i].size())
					requirements.emplace(i, processes[i].size());
		} else {
			// star-root process has not finished yet
			// collect info of what needs to HB the next star-root step
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
			reqNodes.insert(processes[rq.first][ current[rq.first] ]);
		}
		while (!reqNodes.empty()) {
      // find (an arbitrary) head of reqNodes
			auto it = reqNodes.begin();
			assert(it != reqNodes.end());
			const Node *headnd = *it;
			++it;
			while (it != reqNodes.end()) {
        if (hasEdge(*it, headnd, po))
					headnd = *it;
			  ++it;
			}
			// perform the head, replace in reqNodes with its
			// thread successor (if that one is also required)
			result.push_back(headnd->getEvent()->blank_copy());
			unsigned tid = headnd->getProcessID();
			++current[tid];
			assert(reqNodes.count(headnd));
			reqNodes.erase(headnd);
			if (current[tid] < requirements[tid])
				reqNodes.insert(processes[tid][ current[tid] ]);
		}

		if (current[starRoot()] < processes[starRoot()].size()) {
      // perform one step of the star-root process
			result.push_back(processes[starRoot()][ current[starRoot()] ]
											 ->getEvent()->blank_copy());
			++current[starRoot()];
		} else
			done = true;
		
	}
	
  return result;
}
