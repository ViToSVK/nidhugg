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
#include <unordered_set>

#include "Debug.h"
#include "VCBasis.h"
#include "VCAnnotationNeg.h"

using ThreadPairsVclocks = std::vector<std::vector<  std::vector<int>  >>;

using PartialOrder = std::pair<std::unique_ptr<ThreadPairsVclocks>,
                               std::unique_ptr<ThreadPairsVclocks>>;

class VCGraphVclock : public VCBasis {

  // succ[i][j][a] = b means:
  // t_i[<=a] *HB* t_j[b<=]
  // so succ[i][j] fixes a (index of t_i),
  // asks for 'smallest' b (smallest successor within t_j)

  // pred[j][i][b] = a means:
  // t_i[<=a] *HB* t_j[b<=]
  // so pred[j][i] fixes b (index of t_j),
  // asks for 'biggest' a (biggest predecessor within t_i)

  // t_i[a] *HB* t_j[b]
  // bigger a => 'stronger' edge
  // smaller b => 'stronger' edge

  const Node *initial_node;

  std::set<const Node *> nodes;

  std::unordered_map<SymAddrSize, std::unordered_set<const Node *>>
    readsNonroot;

  std::unordered_map<SymAddrSize, std::unordered_set<const Node *>>
    readsRoot;

  std::unordered_map<SymAddrSize, std::unordered_set<const Node *>>
    wNonrootUnord;

  std::unordered_map<SymAddrSize, std::vector<const Node *>> wRoot;

  std::unordered_set<unsigned> leafThreadsWithRorW;

  // [ml][tid][evid] returns idx of first event of thread-tid writing to ml
  // starting from AND INCLUDING evid and going back - (evid, evid-1, .., 0)
  // returns -1 if there is no such write
  std::unordered_map<SymAddrSize, std::vector<std::vector<int>>>
    tw_candidate;

  PartialOrder original;

  // Partial orders that are refinements of original
  std::list<PartialOrder> worklist_ready, worklist_done;

 public:

  const PartialOrder& getOriginal() { return original; }

  bool lessThanTwoLeavesWithRorW() const { return leafThreadsWithRorW.size() < 2; }

  bool empty() const {
    return (processes.empty() && cpid_to_processid.empty() &&
            event_to_node.empty() &&
            !initial_node && nodes.empty() &&
            !original.first.get() && !original.second.get() &&
            readsNonroot.empty() && readsRoot.empty() &&
            wNonrootUnord.empty() && wRoot.empty() &&
            leafThreadsWithRorW.empty() && tw_candidate.empty() &&
            worklist_ready.empty() && worklist_done.empty());
  }

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  ~VCGraphVclock() {
    delete initial_node;
    for (const Node *nd: nodes) {
      delete nd;
    }
  }

  VCGraphVclock() = delete;

  VCGraphVclock(const std::vector<VCEvent>& trace,
                int star_root_index)
  : VCBasis(star_root_index),
    initial_node(new Node(INT_MAX, INT_MAX, nullptr)),
    nodes(),
    readsNonroot(),
    readsRoot(),
    wNonrootUnord(),
    wRoot(),
    leafThreadsWithRorW(),
    tw_candidate(),
    original(std::unique_ptr<ThreadPairsVclocks>(new ThreadPairsVclocks()),
             std::unique_ptr<ThreadPairsVclocks>(new ThreadPairsVclocks())),
    worklist_ready(),
    worklist_done()
      {
        extendGraph(trace, nullptr);
        assert(starRoot() < processes.size() &&
               "Star root index too big (not enough processes in the initial trace)");
      }

  VCGraphVclock(VCGraphVclock&& oth)
  : VCBasis(std::move(oth)),
    initial_node(oth.initial_node),
    nodes(std::move(oth.nodes)),
    readsNonroot(std::move(oth.readsNonroot)),
    readsRoot(std::move(oth.readsRoot)),
    wNonrootUnord(std::move(oth.wNonrootUnord)),
    wRoot(std::move(oth.wRoot)),
    leafThreadsWithRorW(std::move(oth.leafThreadsWithRorW)),
    tw_candidate(std::move(oth.tw_candidate)),
    original(std::move(oth.original)),
    worklist_ready(std::move(oth.worklist_ready)),
    worklist_done(std::move(oth.worklist_done))
      {
        oth.initial_node = nullptr;
        assert(oth.empty());
      }

  VCGraphVclock& operator=(VCGraphVclock&& oth) = delete;

  VCGraphVclock(const VCGraphVclock& oth) = delete;

  // Partial order that will be moved as original
  // Trace and annotation that will extend this copy of the graph
  VCGraphVclock(const VCGraphVclock& oth,
                PartialOrder&& po,
                const std::vector<VCEvent>& trace,
                const VCAnnotation& annotation)
  : VCBasis(oth),
    initial_node(new Node(INT_MAX, INT_MAX, nullptr)),
    nodes(),
    readsNonroot(),
    readsRoot(),
    wNonrootUnord(),
    wRoot(),
    leafThreadsWithRorW(),
    tw_candidate(),
    original(std::move(po)),
    worklist_ready(),
    worklist_done()
      {
        for (const Node *othnd : oth.nodes) {
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
        // lock_vciid_to_node -- fixed below
        // read_vciid_to_node -- fixed below
        extendGraph(trace, &annotation);
      }

  VCGraphVclock& operator=(VCGraphVclock& oth) = delete;

  /* *************************** */
  /* GRAPH EXTENSION             */
  /* *************************** */

 private:

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
  void extendGraph(const std::vector<VCEvent>& trace,
                   const VCAnnotation *annotationPtr);

  /* *************************** */
  /* EDGE QUESTIONS              */
  /* *************************** */

  bool hasEdge(const Node *n1, const Node *n2, const PartialOrder& po) const {
    assert(n1 && n2 && "Do not have such node");
    assert(n1 != n2);
    if (n1 == initial_node)
      return true; // init HB everything
    if (n2 == initial_node)
      return false; // nothing HB init
    if (n1->getProcessID() == n2->getProcessID())
      return n1->getEventID() < n2->getEventID();
    const ThreadPairsVclocks& succ = *(po.first);
    return succ[n1->getProcessID()]
               [n2->getProcessID()]
               [n1->getEventID()] <= (int) n2->getEventID();
  }

  class POcomp {
    const VCGraphVclock& gr;
    const PartialOrder& po;
   public:
    POcomp(const VCGraphVclock& gr, const PartialOrder& po) : gr(gr), po(po) {}
    bool operator() (const Node *n1, const Node *n2) {
      return (gr.hasEdge(n1, n2, po));
    }
  };

  friend class POcomp;

  bool areOrdered(const Node *n1, const Node *n2, const PartialOrder& po) const {
    return hasEdge(n1, n2, po) || hasEdge(n2, n1, po);
  }

  // Given 'nd', returns its first successor within the
  // thread 'thr_id' according to the vector clocks 'succ'
  const Node *getMinSucc(const Node *nd, unsigned thr_id, const PartialOrder& po) const {
    assert(nd && "Do not have such node");
    assert(nd != initial_node && "Asking getMinSucc for the initial node");
    assert(thr_id < processes.size() && "Such thread does not exist");
    assert(nd->getProcessID() != thr_id && "Asking getMinSucc within the same thread");
    const ThreadPairsVclocks& succ = *(po.first);
    int ev_id = succ[nd->getProcessID()]
                    [thr_id]
                    [nd->getEventID()];
    assert(ev_id >= 0);
    return ((unsigned) ev_id < processes[thr_id].size()) ?
      processes[thr_id][ev_id] : nullptr;
  }

  // Given 'nd', returns its latest predecessor within the
  // thread 'thr_id' according to the vector clocks 'pred'
  const Node *getMaxPred(const Node *nd, unsigned thr_id, const PartialOrder& po) const {
    assert(nd && "Do not have such node");
    assert(nd != initial_node && "Asking getMaxPred for the initial node");
    assert(thr_id < processes.size() && "Such thread does not exist");
    assert(nd->getProcessID() != thr_id && "Asking getMaxPred within the same thread");
    const ThreadPairsVclocks& pred = *(po.second);
    int ev_id = pred[nd->getProcessID()]
                    [thr_id]
                    [nd->getEventID()];
    assert(ev_id < (int) processes[thr_id].size());
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
  void addEdge(const Node *n1, const Node *n2, const PartialOrder& po) const;

  // Helper method for addEdge,
  // maintains thread-pair-wise transitivity
  void addEdgeHelp(unsigned ti, unsigned ti_evx,
                   unsigned tj, unsigned tj_evx,
                   const PartialOrder& po) const;

 public:

  /* *************************** */
  /* MAIN ALGORITHM              */
  /* *************************** */

  friend class VCValClosure;

  // Given 'nd', returns a candidate for a tail write from
  // the thread 'thr_id' using 'succ' and 'tw_candidate'
  const Node *getTailWcandidate(const Node *nd, unsigned thr_id, const PartialOrder& po) const;

  // Given 'nd', returns {root tail write, nonroot tail writes} in 'po'
  std::pair<const Node *, std::unordered_set<const Node *>>
    getTailWrites(const Node *nd, const PartialOrder& po) const;

  // Given 'nd', returns a candidates for a head write from the thread 'thr_id' that
  // happens before nd, second are indices from-to for unordered candidate search
  std::pair<const Node *, std::pair<int, int>>
    getHeadWcandidate(const Node *nd, unsigned thr_id, const PartialOrder& po) const;

  // Given 'nd', returns {root head write, nonroot head writes} in 'po'
  std::pair<const Node *, std::unordered_set<const Node *>>
    getHeadWrites(const Node *nd, const PartialOrder& po) const;

  // Used just before ordering non-star-root writes of trace extension
  void initWorklist() {
    assert(worklist_ready.empty());
    assert(worklist_done.empty());
  }

  // Used just before ordering a non-star read-to-be-mutated
  void initWorklist(const PartialOrder& po) {
    assert(worklist_ready.empty());
    assert(worklist_done.empty());
    worklist_done.emplace_back(std::unique_ptr<ThreadPairsVclocks>
                               (new ThreadPairsVclocks(*(po.first))),
                               std::unique_ptr<ThreadPairsVclocks>
                               (new ThreadPairsVclocks(*(po.second))));
  }

  // Used after ordering extension writes or a read-to-be-mutated
  // Each PO will be a target for mutations
  std::list<PartialOrder> dumpDoneWorklist() {
    assert(worklist_ready.empty());

    if (worklist_done.empty())
      return std::list<PartialOrder>();

    assert(!worklist_done.empty());
    auto result = std::list<PartialOrder>();
    result.swap(worklist_done);
    assert(worklist_done.empty());
    assert(!result.empty());

    return result;
  }

  // Given a write node 'nd' and a partial order 'po',
  // determines whether 'nd' is observable in that 'po'
  bool isObservable(const Node *nd, const PartialOrder& po) const;

 private:
  // Helper for isObservable, where the read node is specified
  bool isObservableBy(const Node *writend, const Node *readnd,
                      const PartialOrder& po) const;

 public:

  // Input: partial orders in worklist_ready
  // Output: partial orders in worklist_done
  void orderEventMaz(const VCEvent *ev1, const VCAnnotation& annotation,
                     bool newlyEverGoodWrite, const PartialOrder& po);

  // Returns last nodes of processes that are reads or locks
  std::unordered_set<const Node *> getNodesToMutate(const VCAnnotation& annotation) const {
    auto result = std::unordered_set<const Node *>();
    for (unsigned tid = 0; tid < processes.size(); ++tid) {
      const Node *nd = processes[tid][ processes[tid].size() - 1 ];
      if ((isRead(nd) && !annotation.defines(nd)) ||
          (isLock(nd) && !annotation.isLastLock(nd)))
        result.insert(nd);
    }
    return result;
  }

  // Returns lengths of processes
  std::vector<unsigned> getProcessLengths() const {
    auto result = std::vector<unsigned>();
    result.reserve(processes.size());
    for (unsigned tid = 0; tid < processes.size(); ++tid) {
      assert(processes[tid].size() > 0);
      result.push_back(processes[tid].size() - 1);
    }
    return result;
  }

  // Returns mutation candidates for a read node
  std::map<std::pair<int, VCAnnotation::Loc>, VCAnnotation::Ann>
    getMutationCandidates(const PartialOrder& po,
                          const VCAnnotationNeg& negative, const Node *readnd) const;

  // Linearizes a partial order
  std::vector<VCEvent> linearize(const PartialOrder& po,
                                 const VCAnnotation& annotation) const;

  /* *************************** */
  /* DUMPS                       */
  /* *************************** */

  void dump_po(const PartialOrder& po) const;
  void dump_po() const { dump_po(original); }

  void to_dot(const PartialOrder& po, const char *edge_params = nullptr) const;
  void to_dot(const char *edge_params = nullptr) const { to_dot(original, edge_params); }

};

#endif // _VC_GRAPHVCLOCK_H
