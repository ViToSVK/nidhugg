/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2019 Viktor Toman
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

#ifndef __Z_GRAPH_H__
#define __Z_GRAPH_H__

#include <list>

#include "ZPartialOrder.h"
#include "ZAnnotationNeg.h"


class ZGraph {
 public:
  class Cache {
  public:
    // ML -> thr -> thread-ordered memory-writes
    std::unordered_map
      <SymAddrSize, std::unordered_map
      <unsigned, std::vector<const ZEvent *>>> wm;

    // Read -> its local buffer-write
    std::unordered_map
      <const ZEvent *, const ZEvent *> readWB;

    bool empty() const {
      return (wm.empty() && readWB.empty());
    }
  };

 private:
  ZBasis basis;
  ZPartialOrder po;
  Cache cache;

 public:
  const ZBasis& getBasis() const { return basis; }
  const ZPartialOrder& getPo() const { return po; }
  const Cache& getCache() const { return cache; }

  bool empty() const {
    return (basis.empty() && po.empty() && cache.empty());
  }

  ZPartialOrder copyPO() const {
    return ZPartialOrder(po, basis);
  }


  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  ~ZGraph() {};
  // Empty
  ZGraph();
  // Initial
  ZGraph(const std::vector<ZEvent>& trace, int star_root_index);
  // Moving
  ZGraph(ZGraph&& oth);

  // Partial order that will be moved as original
  // Trace and annotation that will extend this copy of the graph
  ZGraph(const ZGraph& oth,
         ZPartialOrder&& po,
         const std::vector<ZEvent>& trace,
         const ZAnnotation& annotation);

  ZGraph& operator=(ZGraph&& oth) = delete;
  ZGraph(const ZGraph& oth) = delete;
  ZGraph& operator=(const ZGraph& oth) = delete;

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
  //    replace their corresponding 'orig_trace' pointers
  // 2) new events of 'trace' extend existing threads / create new threads
  // 3) partial order is extended to accomodate new
  //    threads+events while keeping all the original info
  // Special case: initial trace extends an empty graph
  void traceToPO(const std::vector<ZEvent>& trace, const ZAnnotation *annotationPtr);


  /* *************************** */
  /* MAIN ALGORITHM              */
  /* *************************** */

 public:
  // In thread thr (and specific aux), starting from ev and going back (ev,ev-1,...,1,0),
  // return first memory-write to ml (resp. index in cache.wm)
  int getTailWindex(const SymAddrSize& ml, unsigned thr, unsigned ev) const;
  const ZEvent *getTailW(const SymAddrSize& ml, unsigned thr, unsigned ev) const;

  // Given 'nd', returns {root tail write, nonroot tail writes} in 'po'
  //std::pair<const Node *, std::unordered_set<const Node *>>
  //getTailWrites(const Node *nd, const PartialOrder& po) const;

  // Given 'nd', returns a candidates for a head write from the thread 'thr_id' that
  // happens before nd, second are indices from-to for unordered candidate search
  //std::pair<const Node *, std::pair<int, int>>
  //getHeadWcandidate(const Node *nd, unsigned thr_id, const PartialOrder& po) const;

  // Given 'nd', returns {root head write, nonroot head writes} in 'po'
  //std::pair<const Node *, std::unordered_set<const Node *>>
  //getHeadWrites(const Node *nd, const PartialOrder& po) const;

  // Given a write node 'nd' and a partial order 'po',
  // determines whether 'nd' is observable in that 'po'
  //bool isObservable(const Node *nd, const PartialOrder& po) const;

 private:
  // Helper for isObservable, where the read node is specified
  //bool isObservableBy(const Node *writend, const Node *readnd,
  //const PartialOrder& po) const;

 public:

  // Input: partial orders in worklist_ready
  // Output: partial orders in worklist_done
  //void orderEventMaz(const ZEvent *ev1, const ZAnnotation& annotation,
  //bool newlyEverGoodWrite, const PartialOrder& po);

  // Returns events-to-mutate in a specified order
  std::list<const ZEvent *> getEventsToMutate(const ZAnnotation& annotation) const;

  // Returns observation candidates for a read node
  std::list<ZObs> getObsCandidates(const ZEvent *read,
                                   const ZAnnotationNeg& negative) const;

  // Linearizes a partial order
  //std::vector<ZEvent> linearize(const PartialOrder& po,
  //const ZAnnotation& annotation) const;

  void dump() const { po.dump(); }

};

#endif // __Z_GRAPH_H__
