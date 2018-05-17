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

VCGraphVclock::~VCGraphVclock()
{
  //for (auto& it : *this)
	//delete &it;
	// TODO

	delete initial_node;
}

VCGraphVclock::VCGraphVclock(const std::vector<VCEvent>& trace)
	: initial_node(new Node(INT_MAX, INT_MAX, nullptr)) {
  /*

  // TraceT is a sequence of VCEvents.
  // (Whatever the TraceT is, the iterators must
  // return the VCEvent& on dereference.)
  template <typename TraceT>
  VCBasis(TraceT& trace, const VCEvent *till = nullptr)
  {
    // mapping from already discovered threads to our processes indexes
    std::unordered_map<IPid, unsigned> mapping;
    mapping.reserve(trace.size() / 2);

    for (const VCEvent& event : trace) {
      IPid pid = event.iid.get_pid();
      unsigned idx;

      auto it = mapping.find(pid);
      if (it == mapping.end()) {
        idx = processes.size();
        mapping[pid] = idx;
        processes.emplace_back(std::vector<const VCEvent *>());
      } else
        idx = it->second;

      processes[idx].push_back(&event);

      // we want processes only till this event?
      // (we know that anything that occured after this event
      //  is irrelevant for us)
      if (till == &event)
        break;
    }
  }


  void reserveMemory() {
    assert(basis.size() > 0);

    // this should be more efficient than to
    // dynamically adjust the hash table
    unsigned reserve_size = 0;
    for (auto&b : basis)
      reserve_size += b.size();
    instr_to_nodeset.reserve(reserve_size);
    mapping.reserve(reserve_size);
  }

  void createInitialNode() {
    initial_node = new Node(basis, INT_MAX, INT_MAX);
    // the initial node is before all other events
    for (int& x : initial_node->successors)
      x = 0;
    mapping.emplace(nullptr, initial_node);
  }

  // create nodes and starting points
  assert(!basis.empty());

  // set of nodes that call pthread_create and pthread_join
  std::vector<Node *> spawns;
  std::vector<Node *> joins;
  spawns.reserve(4);
  joins.reserve(4);

  unsigned process_idx = 0;
  for (auto& base : basis) {
    Node *node;
    Node *last = nullptr;

    // create a new process in nodes
    //  - nodes are basically a copy of basis,
    //  in this moment, but with nodes instead of DCEvents
    processes.emplace_back();

    // create the nodes and add the program structure
    unsigned event_idx = 0;
    for (const DCEvent *event : base) {
      assert(event->iid.get_pid() >= 0);
      node = new Node(basis, process_idx, event_idx);
      mapping.emplace(event, node);
      processes.back().push_back(node);

      if (event->instruction)
        instr_to_node[event->instruction].push_back(node);

      if (last)
        last->addSuccessor(node);

      last = node;

      // XXX: what about function pointer calls?
      if (event->instruction) {
        if (is_function_call(event->instruction, "pthread_create"))
          spawns.push_back(node);
        else if (is_function_call(event->instruction, "pthread_join"))
          joins.push_back(node);
      }

      ++event_idx;
    }

    ++process_idx;
  }

  // fill in the annotation edges
  for (auto& rw_pair : annotation) {
    Node *read = mapAnnotToNode(rw_pair.first);
    Node *write = mapAnnotToNode(rw_pair.second);

    assert(read);
    if (!write)
      write = initial_node;

    write->addSuccessor(read);
  }

  // fill in the threads creation happens-before relation
  for (Node *spwn : spawns) {
    // find the first event of the thread in basis
    for (auto& base : basis) {
      assert(!base.empty());
      if (base[0]->cpid == spwn->getEvent()->childs_cpid) {
        spwn->addSuccessor(getNode(base[0]));
        break;
      }
    }
  }

  for (Node *jn : joins) {
    // find the first event of the thread in basis
    for (auto& base : basis) {
      assert(!base.empty());
      if (base[0]->cpid == jn->getEvent()->childs_cpid) {
        getNode(base[base.size() - 1])->addSuccessor(jn);
        break;
      }
    }
  }

	 */
}

/*
HappensAfterGraphVclock::HappensAfterGraphVclock(HappensAfterGraphVclock&& oth)
: annotation(std::move(oth.annotation)),
  happens_before(std::move(oth.happens_before)),
  basis(oth.basis),
  blocked_before_sat(oth.blocked_before_sat),
  mapping(std::move(oth.mapping)),
  instr_to_node(std::move(oth.instr_to_node)),
  initial_node(oth.initial_node)
{
}

HappensAfterGraphVclock::HappensAfterGraphVclock(const HappensAfterGraphVclock& oth)
: annotation(oth.annotation),
  happens_before(oth.happens_before),
  basis(oth.basis),
  blocked_before_sat(oth.blocked_before_sat)
{
  mapping.reserve(oth.size());
  instr_to_node.reserve(instr_to_node.size());

  for (auto& it : oth.mapping) {
    // copy the node (successors are indices into basis,
    // so the successors will be mapped fine)
    assert(it.second);
    Node *new_nd = new Node(*it.second);

    if (it.first && it.first->instruction)
      instr_to_node[it.first->instruction].emplace_back(new_nd);

    mapping.emplace(it.first, new_nd);
  }

  // also set the initial_node
  initial_node = mapping[nullptr];
  assert(initial_node);
  assert(*this == oth);
}

HappensAfterGraphVclock& HappensAfterGraphVclock::operator=(HappensAfterGraphVclock&& oth) {
  if (&oth != this)
    *this = std::move(oth);
  return *this;
}
*/

bool VCGraphVclock::operator==(const VCGraphVclock& oth) const {
	/*
	if (mapping.size() != oth.mapping.size())
    return false;

  for (auto& it : *this) {
    const Node *oth_nd = oth.getNode(it.getEvent());
    if (!oth_nd  || (it.successors != oth_nd->successors))
      return false;
  }
	*/
  return true; // TODO
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
