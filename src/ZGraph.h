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

    // Leaf ThrID -> memory-writes to chronologically order
    // Leaf ThrID -> annotated reads to chronologically order
    std::unordered_map
      <unsigned, std::list<const ZEvent *>> chrono, chronoAnnR;


    bool empty() const {
      return (wm.empty() && readWB.empty() && chrono.empty());
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
    return ZPartialOrder(po);
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
  // return first memory-write to ml (resp. return its index in cache.wm)
  int getTailWindex(const SymAddrSize& ml, unsigned thr, unsigned ev) const;
  const ZEvent *getTailW(const SymAddrSize& ml, unsigned thr, unsigned ev) const;

  // In thread thr (and specific aux), from all memory-writes
  // conflicting with read that do not happen after the read in partial,
  // return the latest one (resp. return its index in cache.wm)
  int getLatestNotAfterIndex(const ZEvent *read, unsigned thr, const ZPartialOrder& partial) const;
  const ZEvent *getLatestNotAfter(const ZEvent *read, unsigned thr, const ZPartialOrder& partial) const;

  // Get the local buffer-write of the read
  // Returns nullptr if there is no local buffer-write before read
  const ZEvent *getLocalBufferW(const ZEvent *read) const;

  // Returns events-to-mutate in a specified order
  std::list<const ZEvent *> getEventsToMutate(const ZAnnotation& annotation) const;

  // Returns observation candidates for a read node
  std::list<ZObs> getObsCandidates(const ZEvent *read,
                                   const ZAnnotationNeg& negative) const;

  // Returns pairs of po-unordered leaf conflicting events
  // that have to be ordered to get a leaves-chronological-po
  std::vector<std::pair<const ZEvent *, const ZEvent *>>
    chronoOrderPairs(const ZEvent *leafread) const;

  // Linearizes a partial order
  std::vector<ZEvent> linearizeTSO
    (const ZPartialOrder& po, const ZAnnotation& annotation) const;

  void dump() const { po.dump(); }

};

#endif // __Z_GRAPH_H__
