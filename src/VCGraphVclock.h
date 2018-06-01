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

#ifndef _VC_GRAPHVCLOCK_H
#define _VC_GRAPHVCLOCK_H

#include <list>
#include <vector>

#include "VCBasis.h"
#include "VCAnnotation.h"

class VCGraphVclock : public VCBasis {

	// succ[i][j] = Vector clock for t_i -> t_j
	// succ[i][j][a] = b means:
	// t_i[<=a] *HB* t_j[b<=]
	// so succ[i][j] fixes a (index of t_i),
	// asks for 'smallest' b (smallest successor within t_j)
	
	// pred[j][i] = Vector clock for t_j -> t_i
	// pred[j][i][b] = a means:
	// t_i[<=a] *HB* t_j[b<=]
	// so pred[j][i] fixes b (index of t_j),
	// asks for 'biggest' a (biggest predecessor within t_i)

  // t_i[a] *HB* t_j[b]
	// bigger a => 'stronger' edge
	// smaller b => 'stronger' edge
	
  const Node *initial_node;
	std::set<Node *> nodes;

  using ThreadPairsVclocks = std::vector<std::vector<  std::vector<int>  >>;
	
	ThreadPairsVclocks succ_original;
	ThreadPairsVclocks pred_original;

	std::unordered_map<SymAddrSize, std::vector<std::vector<int>>> tw_candidate;

	// Container for pairs succ-pred, each item in the container
	// has some additional partial order upon writes of interest
	std::list< std::pair<ThreadPairsVclocks *, ThreadPairsVclocks *> >
		refinements, // in the value-closing process
		ready_to_mutate, // value-closed
		mutation_refinements; // in the value-closing process after introducing a mutation

	unsigned extension_from;

 public:

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */
	
	~VCGraphVclock() {
		delete initial_node;
		for (Node *nd: nodes) {
      delete nd;
		}
		assert(refinements.empty());
		assert(ready_to_mutate.empty());
		assert(mutation_refinements.empty());
	}
	
	VCGraphVclock() = delete;
	
  VCGraphVclock(const std::vector<VCEvent>& trace)
	: VCBasis(),
		initial_node(new Node(INT_MAX, INT_MAX, nullptr))
			{ extendGraph(trace); }

  VCGraphVclock(VCGraphVclock&& oth) = default;

  VCGraphVclock& operator=(VCGraphVclock&& oth) = delete;

  VCGraphVclock(const VCGraphVclock& oth) = delete;

	// Pointers to Vclocks that will be copied as succ/pred_original
	// Trace that will extend this copy of the graph
  VCGraphVclock(const VCGraphVclock& oth,
							  const ThreadPairsVclocks *succ,
							  const ThreadPairsVclocks *pred,
								const std::vector<VCEvent>& trace)
	: VCBasis(oth),
		initial_node(new Node(INT_MAX, INT_MAX, nullptr)),
    succ_original(*succ),
    pred_original(*pred),
		tw_candidate(),
		refinements(),
		ready_to_mutate(),
		mutation_refinements()
			{
				for (Node *othnd : oth.nodes) {
					Node *nd = new Node(*othnd);
					nodes.insert(nd);
					assert(nd->getProcessID() < processes.size() &&
								 nd->getEventID() < processes[nd->getProcessID()].size() &&
								 processes[nd->getProcessID()][nd->getEventID()] == othnd);
					processes[nd->getProcessID()][nd->getEventID()] = nd;
				}
				// nodes -- fixed above
				// processes -- fixed above
				// cpid_to_processid -- no need to fix
				// event_to_node -- fixed below
				// read_vciid_to_node -- fixed below
				extendGraph(trace);
	    }
	
  VCGraphVclock& operator=(VCGraphVclock& oth) = delete;

  /* *************************** */
  /* GRAPH EXTENSION             */
  /* *************************** */

	// At the point of calling the method, this graph
	// is linked to some 'orig_trace', graph's nodes
	// point to the events of 'orig_trace'
	// Argument: 'trace' - extension of 'orig_trace'
	// The method links this graph with 'trace' as follows:
	// 1) events of 'trace' that already happened in 'orig_trace'
	//    are linked to the corresponding (already existing) nodes
	// 2) new events of 'trace' create new nodes and in the basis
	//    they extend existing threads / create new threads
	// 3) succ/pred_original are extended to accomodate new
	//    threads+nodes while keeping all the original info
	// Special case: initial trace extends an empty graph
	void extendGraph(const std::vector<VCEvent>& trace);
	
  /* *************************** */
  /* EDGE QUESTIONS              */
  /* *************************** */

  bool hasEdge(const Node *n1, const Node *n2,
							 const ThreadPairsVclocks *succ) const {
    assert(n1 && n2 && "Do not have such node");
		assert(n1 != n2);
		if (n1 == initial_node)
			return true; // init HB everything
    if (n2 == initial_node)
      return false; // nothing HB init
		if (n1->getProcessID() == n2->getProcessID())
			return n1->getEventID() < n2->getEventID();
    return (*succ)[n1->getProcessID()]
		              [n2->getProcessID()]
			            [n1->getEventID()] <= (int) n2->getEventID();
  }

  bool areOrdered(const Node *n1, const Node *n2,
									const ThreadPairsVclocks *succ) const {
    return hasEdge(n1, n2, succ) || hasEdge(n2, n1, succ);
  }

	// Given 'nd', returns its first successor within the
	// thread 'thr_id' according to the vector clocks 'succ'
	const Node *getMinSucc(const Node *nd, unsigned thr_id,
						     const ThreadPairsVclocks *succ) const {
		assert(nd && "Do not have such node");
		assert(nd != initial_node && "Asking getMinSucc for the initial node");
		assert(thr_id < processes.size() && "Such thread does not exist");
		assert(nd->getProcessID() != thr_id && "Asking getMinSucc within the same thread");
		int ev_id = (*succ)[nd->getProcessID()]
		                   [thr_id]
		                   [nd->getEventID()];
		assert(ev_id >= 0);
		return ((unsigned) ev_id < processes[thr_id].size()) ?
			processes[thr_id][ev_id] : nullptr;
	}

	// Given 'nd', returns its latest predecessor within the
	// thread 'thr_id' according to the vector clocks 'pred'
	const Node *getMaxPred(const Node *nd, unsigned thr_id,
								 const ThreadPairsVclocks *pred) const {
		assert(nd && "Do not have such node");
		assert(nd != initial_node && "Asking getMaxPred for the initial node");
		assert(thr_id < processes.size() && "Such thread does not exist");
		assert(nd->getProcessID() != thr_id && "Asking getMaxPred within the same thread");
		int ev_id = (*pred)[nd->getProcessID()]
			                 [thr_id]
			                 [nd->getEventID()];
		assert(ev_id < processes[thr_id].size());
		return (ev_id >= 0) ?
			processes[thr_id][ev_id] : nullptr;
	}

  /* *************************** */
  /* EDGE ADDITION               */
  /* *************************** */

	// Adds an edge between two unordered nodes
	// This method maintains:
	// 1) thread-pair-wise transitivity
	// 2) complete transitivity
  void addEdge(const Node *n1, const Node *n2,
							 ThreadPairsVclocks *succ,
							 ThreadPairsVclocks *pred);

 private:

	// Helper method for addEdge,
	// maintains thread-pair-wise transitivity
	void addEdgeHelp(unsigned ti, unsigned ti_evx,
									 unsigned tj, unsigned tj_evx,
									 ThreadPairsVclocks *succ,
									 ThreadPairsVclocks *pred);

 public:
	
  // TODO: alg. helpers

  void to_dot(const char *edge_params=nullptr) const;
	
};

#endif // _VC_GRAPHVCLOCK_H
