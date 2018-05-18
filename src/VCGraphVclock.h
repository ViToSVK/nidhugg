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

#include <vector>

#include "VCBasis.h"

class VCGraphVclock : public VCBasis {

  const Node *initial_node;
	std::set<Node *> nodes;
	
	// edges[i][j] = Vector clock for t_i -> t_j
	// edges[i][j][k] = x means:
	// t_i[<=k] *HB* t_j[x<=]
	std::vector<std::vector<  std::vector<unsigned>  >> edges;

 public:
	~VCGraphVclock() {
		delete initial_node;
		for (Node *nd: nodes)
			delete nd;
	}
	
	VCGraphVclock() = delete;
  VCGraphVclock(const std::vector<VCEvent>& trace);

  VCGraphVclock(VCGraphVclock&& oth)
	: VCBasis(std::move(oth)),
		initial_node(std::move(oth.initial_node)),
		nodes(std::move(oth.nodes)),
		edges(std::move(oth.edges)) {}
		
  VCGraphVclock& operator=(VCGraphVclock&& oth) {
		if (&oth != this) {
			delete initial_node;
			for (Node *nd: nodes)
				delete nd;
			initial_node = std::move(oth.initial_node);
			nodes = std::move(oth.nodes);
			edges = std::move(oth.edges);
			VCBasis::operator=(std::move(oth));
		}
		return *this;
	};

  VCGraphVclock(const VCGraphVclock& oth)
	: VCBasis(oth),
		initial_node(new Node(INT_MAX, INT_MAX, nullptr)),
		edges(oth.edges) {
		for (Node *nd : oth.nodes)
			nodes.insert(new Node(*nd));
	}
	
  VCGraphVclock& operator=(VCGraphVclock& oth) = delete;

  bool hasEdge(const VCEvent& e1, const VCEvent& e2) const {
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    return hasEdge(n1, n2);
  }

  bool hasEdge(const Node *n1, const Node *n2) const {
    assert(n1 && n2 && "Do not have such node");
		if (n1 == initial_node)
			return true; // init HB everything
    if (n2 == initial_node)
      return false; // nothing HB init
		if (n1->getProcessID() == n2->getProcessID())
			return n1->getEventID() < n2->getEventID();
    return edges[n1->getProcessID()]
		            [n2->getProcessID()]
		            [n1->getEventID()] <= n2->getEventID();
  }

  bool areOrdered(const VCEvent& e1, const VCEvent& e2) const {
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    return areOrdered(n1, n2);
  }

  bool areOrdered(const Node *n1, const Node *n2) const {
    return hasEdge(n1, n2) || hasEdge(n2, n1);
  }

  void addEdge(const VCEvent& e1, const VCEvent& e2) {
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    addEdge(n1, n2);
  }
	
  void addEdge(const Node *n1, const Node *n2) {
    assert(n1 && n2 && "Do not have such node");
		assert(n1 != initial_node && n2 != initial_node
					 && "Can not add an edge from/to the initial node");
		int i = n1->getProcessID();
		int j = n2->getProcessID();
		assert(i != j && "Can not add an edge within the same process");
		unsigned x = n2->getEventID();
		for (int k=n1->getEventID(); k>=0; --k) {
      if (edges[i][j][k] <= x)
				break; // since for smaller k's also <= x
			edges[i][j][k] = x;
		}
  }

  // Returns {true, order} if the graph is acyclic and
  // {false, _} if there's a cycle
  //std::pair<bool, std::vector<Node *>> computeTopoOrder();
	
  // Linearize the PO induced by this graph to create a trace
  //std::vector<VCEvent> linearize();

  void to_dot(const char *edge_params=nullptr) const;
  void dump() const;

};

#endif // _VC_GRAPHVCLOCK_H
