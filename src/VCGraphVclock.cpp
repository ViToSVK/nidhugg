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

#include "VCTBHelpers.h" // isWrite
#include "VCGraphVclock.h"
#include "SCC.h"

/* *************************** */
/* CONSTRUCTORS                */
/* *************************** */

VCGraphVclock::VCGraphVclock(const std::vector<VCEvent>& trace)
	: VCBasis(), initial_node(new Node(INT_MAX, INT_MAX, nullptr)) {

  assert(!trace.empty());
	
	// Faster than cpid_to_processid
	std::unordered_map<int, unsigned> ipid_mapping;
  ipid_mapping.reserve(trace.size() / 2);
	cpid_to_processid.reserve(trace.size() / 2);
	processes.reserve(trace.size() / 2);
	event_to_node.reserve(trace.size());
	read_vciid_to_node.reserve(trace.size() / 2);

  std::vector<Node *> spawns;
  std::vector<Node *> joins;
  spawns.reserve(8);
  joins.reserve(8);  

  for (auto traceit = trace.begin(); traceit != trace.end(); ++traceit) {
		const VCEvent *ev = &(*traceit);
    int pid = ev->iid.get_pid();
    
    unsigned proc_idx, ev_idx;
		
		auto ipidit = ipid_mapping.find(pid);
		if (ipidit == ipid_mapping.end()) {
			// First event of a new process
			proc_idx = processes.size();
			ev_idx = 0;
			processes.emplace_back(std::vector<Node *>());
			ipid_mapping.emplace(pid, proc_idx);
			cpid_to_processid.emplace(ev->cpid, proc_idx);
		} else {
      // Another event of a known process
			proc_idx = ipidit->second;
			ev_idx = processes[proc_idx].size();
		}
		
    Node *nd = new Node(proc_idx, ev_idx, ev);

		nodes.insert(nd);
		processes[proc_idx].push_back(nd);
		assert(processes[proc_idx].size() == ev_idx + 1);
		event_to_node.emplace(ev, nd);
		if (isRead(*ev))
			read_vciid_to_node.emplace(VCIID(*ev), nd);

		// XXX: what about function pointer calls?
		if (ev->instruction) {
			if (is_function_call(ev->instruction, "pthread_create"))
				spawns.push_back(nd);
			else if (is_function_call(ev->instruction, "pthread_join"))
				joins.push_back(nd);
	  }
		
	}
	
	// EDGES - initialize
	succ_original.reserve(processes.size());
	pred_original.reserve(processes.size());
	for (unsigned i=0; i<processes.size(); i++) {
    succ_original.push_back(std::vector<std::vector<int>>());
		succ_original[i].reserve(processes.size());
    pred_original.push_back(std::vector<std::vector<int>>());
		pred_original[i].reserve(processes.size());

		for (unsigned j=0; j<processes.size(); j++) {
      succ_original[i].push_back(std::vector<int>());
			succ_original[i][j].reserve(i == j ? 0 : processes[i].size());
      pred_original[i].push_back(std::vector<int>());
			pred_original[i][j].reserve(i == j ? 0 : processes[i].size());

			if (i != j)
				for (unsigned k=0; k<processes[i].size(); k++) {
					succ_original[i][j][k] = INT_MAX;
					pred_original[i][j][k] = -1;
				}
		}
	}

	// EDGES - spawns
	for (Node *spwn : spawns) {
		int spwn_pid = cpidToProcessID(spwn->getEvent()->childs_cpid);
		assert(spwn_pid >= 0 && spwn_pid < processes.size());
    bool res = addEdge(spwn, processes[spwn_pid][0],
											 &succ_original, &pred_original);
		assert(!res && "Adding spawn edges created a cycle");
		((void)(res)); // so res does not appear unused on release
	}

	// EDGES - joins
	for (Node *jn : joins) {
		int jn_pid = cpidToProcessID(jn->getEvent()->childs_cpid);
		assert(jn_pid >= 0 && jn_pid < processes.size());
    bool res = addEdge(processes[jn_pid][processes[jn_pid].size() - 1], jn,
											 &succ_original, &pred_original);
		assert(!res && "Adding join edges created a cycle");
		((void)(res)); // so res does not appear unused on release		
	}

}

/* *************************** */
/* GRAPH EXTENSION             */
/* *************************** */

// Extends this graph so it corresponds to 'trace'
// Checks the header file for the method description
void VCGraphVclock::extendGraph(const std::vector<VCEvent>& trace) {
  // TBD
}

/* *************************** */
/* EDGE ADDITION               */
/* *************************** */

// Returns true iff adding the edge would create a cycle
// This method maintains:
// 1) thread-pair-wise transitivity
// 2) complete transitivity
bool VCGraphVclock::addEdge(const Node *n1, const Node *n2,
														ThreadPairsVclocks *succ,
														ThreadPairsVclocks *pred) {
	assert(n1 && n2 && "Do not have such node");
	assert(n1 != initial_node && n2 != initial_node
				 && "Can not add an edge from/to the initial node");
	unsigned ti = n1->getProcessID();
	unsigned tj = n2->getProcessID();
	assert(ti != tj && "Can not add an edge within the same process");
	unsigned ti_evx = n1->getEventID();
	unsigned tj_evx = n2->getEventID();
	assert(ti_evx < processes[ti].size() && tj_evx < processes[tj].size()
				 && "Invalid node indices");
	
	std::pair<bool, bool> addEdgeResult =
									addEdgeHelp(ti, ti_evx, tj, tj_evx, succ, pred);
	
  if (addEdgeResult.second) {
    assert(!addEdgeResult.first);
		return true;
	}

  // Maintenance of complete transitivity
  // Collect nodes from different threads with edges:
	// 1) to   ti[ti_evx]
	// 2) from tj[tj_evx]
	
	std::set< std::pair<unsigned, unsigned>> before_tiEvent, after_tjEvent;
	for (unsigned k = 0; k<processes.size(); ++k) {
    if (k != ti && k != tj) {
      int maxbefore = (*pred)[ti][k][ti_evx];
			if (maxbefore >= 0)
				before_tiEvent.insert(std::pair<unsigned, unsigned>
															(k, maxbefore));

			int minafter = (*succ)[tj][k][tj_evx];
			assert(minafter >= 0);
			if ((unsigned) minafter < processes[k].size())
				after_tjEvent.insert(std::pair<unsigned, unsigned>
														 (k, minafter));
		}
	}

	// Try to add an edge between
	// each of 1) and each of 2)
	// (if they belong to different threads)

	for (auto& bef : before_tiEvent)
		for (auto& aft : after_tjEvent)
			if (bef.first != aft.first) {
        addEdgeResult = addEdgeHelp(bef.first, bef.second,
																		aft.first, aft.second,
																		succ, pred);
				if (addEdgeResult.second) {
					assert(!addEdgeResult.first);
					// careful: even though this edge was not added
					// in order not to create a cycle, the creation
					// of this cycle is already unavoidable in the graph
					// (ie transitively closing the graph has to create
					// this cycle), so the graph is already useless for us
					return true;
				}
			}

	return false;
}

// first)  true iff the edge was added
// second) true iff adding the edge would create a cycle
// This method maintains thread-pair-wise transitivity
std::pair<bool, bool> VCGraphVclock::addEdgeHelp(unsigned ti, unsigned ti_evx,
																								 unsigned tj, unsigned tj_evx,
																								 ThreadPairsVclocks *succ,
																								 ThreadPairsVclocks *pred) {
  if ((*succ)[ti][tj][ti_evx] <= (int) tj_evx) {
    // same/stronger edge already present
		assert((*pred)[tj][ti][tj_evx] >= (int) ti_evx
					 && "Inconsistent succ/pred vector clocks");
		return std::pair<bool, bool>(false, false);
	}

	// weaker edge present, will add this one
	assert((*pred)[tj][ti][tj_evx] < (int) ti_evx
				 && "Inconsistent succ/pred vector clocks");

	// check for cycle
  if ((*succ)[tj][ti][tj_evx] <= (int) ti_evx) {
		assert((*pred)[ti][tj][ti_evx] >= (int) tj_evx
					 && "Inconsistent succ/pred vector clocks");		
    return std::pair<bool, bool>(false, true);
	}

	// adding the edge will not create a cycle
  // checking for cycle during maintenance
	// of thread-pair-wise is unnecessary, since
	// a cycle would have been found above already

	(*succ)[ti][tj][ti_evx] = tj_evx;
	(*pred)[tj][ti][tj_evx] = ti_evx;	
	
  // maintenance of transitivity of succ[ti][tj]
	// tj_evx is fixed, everything in ti before ti_evx
	// also happens-before tj_evx
	
	for (int ti_before_evx = ti_evx - 1;
			 ti_before_evx >= 0;
			 --ti_before_evx) {
		if ((*succ)[ti][tj][ti_before_evx] <= (int) tj_evx)
      break; // since for smaller indices also <= tj_evx
		(*succ)[ti][tj][ti_before_evx] = tj_evx;
	}

  // maintenance of transitivity of pred[tj][ti]
	// ti_evx is fixed, everything in tj after tj_evx
	// also happens-after ti_evx
	
	for (unsigned tj_after_evx = tj_evx + 1;
			 tj_after_evx < processes[tj].size();
			 ++tj_after_evx) {
		if ((*pred)[tj][ti][tj_after_evx] >= (int) ti_evx)
			break; // since for bigger indices also >= ti_evx
		(*pred)[tj][ti][tj_after_evx] = ti_evx;
	}

	return std::pair<bool, bool>(true, false);
}

// Returns true iff performing the closure implies creating a cycle
bool VCGraphVclock::closure(ThreadPairsVclocks *succ,
														ThreadPairsVclocks *pred) {
	// TBD: Add closure rules here
  return false;
}
