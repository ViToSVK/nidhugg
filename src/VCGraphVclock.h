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

class VCGraphVclock : public VCBasis {

  const Node *initial_node;
	std::set<Node *> nodes;

  using ThreadPairsVclocks = std::vector<std::vector<  std::vector<int>  >>;
	
	// succ[i][j] = Vector clock for t_i -> t_j
	// succ[i][j][a] = b means:
	// t_i[<=a] *HB* t_j[b<=]
	// so succ[i][j] fixes a (index of t_i),
	// asks for 'smallest' b (smallest successor within t_j)
	ThreadPairsVclocks succ_original;

	// pred[j][i] = Vector clock for t_j -> t_i
	// pred[j][i][b] = a means:
	// t_i[<=a] *HB* t_j[b<=]
	// so pred[j][i] fixes b (index of t_j),
	// asks for 'biggest' a (biggest predecessor within t_i)
	ThreadPairsVclocks pred_original;

  // t_i[a] *HB* t_j[b]
	// bigger a => 'stronger' edge
	// smaller b => 'stronger' edge

	// Container for pairs succ-pred, each item in the container
	// has some additional partial order upon writes of interest
	std::list< std::pair<ThreadPairsVclocks *, ThreadPairsVclocks *> >
		edges_worklist;

 public:

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */
	
	~VCGraphVclock() {
		delete initial_node;
		for (Node *nd: nodes) {
      delete nd;
		}
		for (auto succpredpair : edges_worklist) {
			delete succpredpair.first;
			delete succpredpair.second;
		}
	}
	
	VCGraphVclock() = delete;
  VCGraphVclock(const std::vector<VCEvent>& trace);

  VCGraphVclock(VCGraphVclock&& oth)
	: VCBasis(std::move(oth)),
		initial_node(std::move(oth.initial_node)),
		nodes(std::move(oth.nodes)),
		succ_original(std::move(oth.succ_original)),
		pred_original(std::move(oth.pred_original)),
		edges_worklist(std::move(oth.edges_worklist)) {}
		
  VCGraphVclock& operator=(VCGraphVclock&& oth) {
		if (&oth != this) {
			// destroy visible resources of lhs
			delete initial_node;
			for (Node *nd: nodes) {
				delete nd;
			}
			for (auto succpredpair : edges_worklist) {
				delete succpredpair.first;
				delete succpredpair.second;
			}
			// move assign bases+members from rhs
			initial_node = std::move(oth.initial_node);
			nodes = std::move(oth.nodes);
			succ_original = std::move(oth.succ_original);
			pred_original = std::move(oth.pred_original);
			edges_worklist = std::move(oth.edges_worklist);
			VCBasis::operator=(std::move(oth));
			// rhs is now resource-less
		}
		return *this;
	};

  VCGraphVclock(const VCGraphVclock& oth) = delete;

	// Do NOT copy edges_worklist, provide pointers to
	// TPVclocks that will be copied as the succ/pred_original
  VCGraphVclock(const VCGraphVclock& oth,
							  const ThreadPairsVclocks *succ,
							  const ThreadPairsVclocks *pred)
	: VCBasis(oth),
		initial_node(new Node(INT_MAX, INT_MAX, nullptr)),
    succ_original(*succ),
    pred_original(*pred),
		edges_worklist() {
		for (Node *nd : oth.nodes)
			nodes.insert(new Node(*nd));
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
	void extendGraph(const std::vector<VCEvent>& trace);
	
  /* *************************** */
  /* EDGE QUESTIONS              */
  /* *************************** */
	
  bool hasEdge(const VCEvent& e1, const VCEvent& e2,
							 const ThreadPairsVclocks *succ) const {
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    return hasEdge(n1, n2, succ);
  }

  bool hasEdge(const Node *n1, const Node *n2,
							 const ThreadPairsVclocks *succ) const {
    assert(n1 && n2 && "Do not have such node");
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

  bool areOrdered(const VCEvent& e1, const VCEvent& e2,
									const ThreadPairsVclocks *succ) const {
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    return areOrdered(n1, n2, succ);
  }

  bool areOrdered(const Node *n1, const Node *n2,
									const ThreadPairsVclocks *succ) const {
    return hasEdge(n1, n2, succ) || hasEdge(n2, n1, succ);
  }

	// Given 'nd', returns its first successor within the
	// thread 'thr_id' according to the vector clocks 'succ'
	int getMinSucc(const Node *nd, unsigned thr_id,
						     const ThreadPairsVclocks *succ) const {
		assert(nd && "Do not have such node");
		assert(nd != initial_node && "Asking getMinSucc for the initial node");
		assert(thr_id < processes.size() && "Such thread does not exist");
		return (*succ)[nd->getProcessID()]
			            [thr_id]
			            [nd->getEventID()];
	}

	// Given 'nd', returns its latest predecessor within the
	// thread 'thr_id' according to the vector clocks 'pred'
	int getMaxPred(const Node *nd, unsigned thr_id,
								 const ThreadPairsVclocks *pred) const {
		assert(nd && "Do not have such node");
		assert(nd != initial_node && "Asking getMaxPred for the initial node");
		assert(thr_id < processes.size() && "Such thread does not exist");
		return (*pred)[nd->getProcessID()]
			            [thr_id]
			            [nd->getEventID()];
	}

  /* *************************** */
  /* EDGE ADDITION               */
  /* *************************** */

	// FUTURE BELOW
  //using edges = std::set<std::pair<const Node *, const Node *>>;
	// <0> 'Strong' added edges (ie maximal wrt edge strength ordering)
	// <1> All added edges (including those hidden by Vclock principle)
	// <2> true iff a cycle appeared
  // FUTURE END

	bool addEdge(const VCEvent& e1, const VCEvent& e2,
							 ThreadPairsVclocks *succ,
							 ThreadPairsVclocks *pred) {
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    return addEdge(n1, n2, succ, pred);
  }

	// Returns true iff adding the edge would create a cycle
	// This method maintains:
	// 1) thread-pair-wise transitivity
	// 2) complete transitivity
  bool addEdge(const Node *n1, const Node *n2,
							 ThreadPairsVclocks *succ,
							 ThreadPairsVclocks *pred);

	// Returns true iff performing the closure implies creating a cycle
  bool closure(ThreadPairsVclocks *succ,
							 ThreadPairsVclocks *pred);
	
  void to_dot(const char *edge_params=nullptr) const;
  void dump() const;

 private:

	// first)  true iff the edge was added
	// second) true iff adding the edge would create a cycle
	// This method maintains thread-pair-wise transitivity
	std::pair<bool, bool> addEdgeHelp(unsigned ti, unsigned ti_evx,
																		unsigned tj, unsigned tj_evx,
																		ThreadPairsVclocks *succ,
																		ThreadPairsVclocks *pred);

};

#endif // _VC_GRAPHVCLOCK_H
