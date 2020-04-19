/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2020 Viktor Toman
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


typedef std::vector<const ZEvent *> LineT;
typedef std::vector<LineT> LinesT;


class ZGraph {
  friend class ZPartialOrder;
  // TSO/PSO
 public:
  const bool tso;


  // LINES (changed at recursion child with new ZEvent pointers from new trace)
 private:
  const ZEvent init; ////
  LinesT lines; ////
 public:
  const ZEvent *initial() const { return &init; }
  const LineT& operator()(std::pair<unsigned, int> ids) const;
  const LineT& operator()(unsigned thread_id, int aux_id) const;
  //
  const ZEvent *getEvent(unsigned thread_id, int aux_id, unsigned event_id) const;
  const ZEvent *getEvent(const ZObs& obs) const;
  const ZEvent *getUnlockOfThisLock(const ZObs& obs) const;
  void addLine(const ZEvent *ev);
  void addEvent(const ZEvent *ev);
  void replaceEvent(const ZEvent *oldEv, const ZEvent *newEv);
  void shrink();


  // THREAD_AUX->LINE_ID (retained accross recursion children)
 private:
  std::unordered_map<std::pair<unsigned, int>, unsigned> thread_aux_to_line_id; ////
  std::vector<std::set<int>> threads_auxes; ////
  std::unordered_map<unsigned, std::vector<SymAddrSize>> pso_thr_mlauxes; ////
  unsigned lineID(unsigned thread_id, int aux_id) const;
  unsigned lineID(const ZEvent *ev) const;
 public:
  bool hasThreadAux(std::pair<unsigned, int> ids) const;
  bool hasThreadAux(unsigned thread_id, int aux_id) const;
  std::vector<unsigned> real_sizes_minus_one() const;
  unsigned number_of_threads() const;
  // All auxiliary indices for thread_id
  const std::set<int>& auxes(unsigned thread_id) const;
  // Auxiliary index for the ml in thread thr
  // Returns -1 if thr has no writes
  // Otherwise, the answer is always 0 in TSO
  int auxForMl(const SymAddrSize& ml, unsigned thr) const;
  int psoGetAux(const ZEvent* writeM);

  // PROCSEQ->THREAD_ID (retained accross recursion children)
 private:
  std::unordered_map<std::vector<int>, unsigned> proc_seq_to_thread_id; ////
 public:
  // <threadID, added_with_this_call?>
  std::pair<unsigned, bool> getThreadID(const std::vector<int>& proc_seq);
  std::pair<unsigned, bool> getThreadID(const ZEvent * ev);
  unsigned getThreadIDnoAdd(const std::vector<int>& proc_seq) const;
  unsigned getThreadIDnoAdd(const ZEvent * ev) const;


  // EVENT->POSITION (changed at recursion child with new ZEvent pointers from new trace)
 private:
  std::unordered_map<const ZEvent *, std::pair<unsigned, unsigned>> event_to_position; ////
 public:
  bool hasEvent(const ZEvent *ev) const;


  // PO
 private:
  ZPartialOrder po;
 public:
  const ZPartialOrder& getPo() const { return po; }
  ZPartialOrder copyPO() const { return ZPartialOrder(po); }


  // CACHE
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
  Cache cache;
 public:
  const Cache& getCache() const { return cache; }


  bool empty() const { return (lines.empty() && po.empty() && cache.empty()); }
  size_t size() const { return lines.size(); }
  size_t events_size() const {
    size_t res = 0;
    for (auto& ln : lines)
      res += ln.size();
    return res;
  }

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  ~ZGraph() {};
  // Empty
  ZGraph();
  // Initial
  ZGraph(const std::vector<ZEvent>& trace, bool tso);
  // Moving
  ZGraph(ZGraph&& oth);

  // Extending
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
  int getTailWindex(const SymAddrSize& ml, unsigned thr, int evX) const;
  const ZEvent *getTailW(const SymAddrSize& ml, unsigned thr, int ev) const;

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


  std::string to_string() const { return po.to_string(); }
  void dump() const { po.dump(); }

};

#endif // __Z_GRAPH_H__
