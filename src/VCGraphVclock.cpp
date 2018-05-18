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

VCGraphVclock::VCGraphVclock(const std::vector<VCEvent>& trace)
	: VCBasis(), initial_node(new Node(INT_MAX, INT_MAX, nullptr)) {

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
	edges.reserve(processes.size());
	for (unsigned i=0; i<processes.size(); i++) {
    edges.push_back(std::vector<std::vector<unsigned>>());
		edges[i].reserve(processes.size());

		for (unsigned j=0; j<processes.size(); j++) {
      edges[i].push_back(std::vector<unsigned>());
			edges[i][j].reserve(i == j ? 0 : processes[i].size());

			for (unsigned k=0; k<processes[i].size(); k++)
				edges[i][j][k] = INT_MAX;
		}
	}

	// EDGES - spawns
	for (Node *spwn : spawns) {
		int spwn_pid = cpidToProcessID(spwn->getEvent()->childs_cpid);
		assert(spwn_pid >= 0 && spwn_pid < processes.size());
    addEdge(spwn, processes[spwn_pid][0]);
	}

	// EDGES - joins
	for (Node *jn : joins) {
		int jn_pid = cpidToProcessID(jn->getEvent()->childs_cpid);
		assert(jn_pid >= 0 && jn_pid < processes.size());
    addEdge(processes[jn_pid][processes[jn_pid].size() - 1], jn);
	}

}



/*
std::pair<bool, std::vector<Node *>>
VCGraphVclock::computeTopoOrder()
{
  struct info {
    unsigned dfsid = 0, lowpt = 0;
    bool is_on_stack = false;
  };


  // initialize the info about nodes
  info init_info; // info of the initial node
  std::vector<std::vector<info>> infos;
  infos.resize(processes.size());
  unsigned i = 0;
  for (auto& inf : infos)
    inf.resize(processes[i++].size());

  std::function<info&(Node *)> get_info = [&infos, &init_info](Node *nd)-> info& {
    if (nd->process_id == INT_MAX)
      return init_info;
    return infos[nd->process_id][nd->event_id];
  };

  std::stack<Node *> stack;
  std::vector<Node *> ret;
  unsigned index = 0;

  std::function<bool(Node *n)> _compute = [&](Node *n) -> bool {
    ++index;
    auto& inf = get_info(n);
    inf.dfsid = index;
    inf.lowpt = index;
    stack.push(n);
    inf.is_on_stack = true;

    for (auto it = succ_begin(n), end = succ_end(n); it != end; ++it) {
      Node *succ = *it;
      auto& succ_inf = get_info(succ);
      if (succ_inf.dfsid == 0) {
        assert(succ_inf.is_on_stack == false);
        if (!_compute(succ))
          return false;
        inf.lowpt = std::min(inf.lowpt, succ_inf.lowpt);
      } else if (succ_inf.is_on_stack) {
        return false; // cycle
      }
    }

    if (inf.lowpt == inf.dfsid) {
        int its = 0;
        while (get_info(stack.top()).dfsid >= inf.dfsid) {
          ++its;
          Node *w = stack.top();
          stack.pop();

          get_info(w).is_on_stack = false;
          ret.push_back(w);

          if (stack.empty())
              break;
        }
        assert(its == 1);
    }

    return true;
  };

  if (!_compute(initial_node))
    return {false, {}};

  assert(stack.empty());

  return {true, ret};
}


// this works since the graph is acyclic
std::vector<VCEvent> VCGraphVclock::linearize()
{
  std::vector<VCEvent> trace;

  auto order = computeTopoOrder();
  assert(order.first && "The graph is not acyclic");

  // the order is reversed
  assert(order.second.size() == size());
  for (unsigned i = order.second.size(); i > 0; --i)
    if (order.second[i - 1]->getEvent()) {
      trace.emplace_back(order.second[i - 1]->getEvent()->blank_copy());
      trace.back().id = trace.size() - 1;
    } else
      // the initial event must be the first event
      assert(i == size());

  return trace;
}
*/
