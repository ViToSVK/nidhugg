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


enum class MemoryModel {
  SC, TSO, PSO
};


typedef std::vector<const ZEvent *> LineT;
typedef std::vector<LineT> LinesT;


class ZGraph {
  friend class ZPartialOrder;
 public:
  const MemoryModel model;

  // LINES (changed at recursion child with new ZEvent pointers from new trace)
 private:
  const ZEvent init; ////
  LinesT lines; ////
 public:
  const ZEvent * const initial() const { return &init; }
  const LineT& operator()(std::pair<unsigned, int> ids) const;
  const LineT& operator()(unsigned thread_id, int aux_id) const;
  //
  const ZEvent * const getEvent(unsigned thread_id, int aux_id, unsigned event_id) const;
  const ZEvent * const getEvent(const ZEventID& id) const;
  void addLine(const ZEvent *ev);
  void addEvent(const ZEvent *ev);
  void shrink();

  // THREAD_AUX->LINE_ID (retained accross recursion children)
 private:
  std::unordered_map<std::pair<unsigned, int>, unsigned> thread_aux_to_line_id; ////
  std::map<int, std::set<int>> threads_auxes; ////
  unsigned lineID(unsigned thread_id, int aux_id) const;
  unsigned lineID(const ZEvent *ev) const;
 public:
  bool hasThreadAux(std::pair<unsigned, int> ids) const;
  bool hasThreadAux(unsigned thread_id, int aux_id) const;
  unsigned number_of_threads() const;
  std::set<unsigned> get_threads() const;
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
  void addThreadID(const ZEvent *ev);
  std::pair<unsigned, bool> getThreadID(const ZEvent * ev);
  unsigned getThreadIDnoAdd(const std::vector<int>& proc_seq) const;
  unsigned getThreadIDnoAdd(const ZEvent * ev) const;
  std::unordered_set<std::vector<int>> all_proc_seq() const;

  // EVENT->POSITION (changed at recursion child with new ZEvent pointers from new trace)
 private:
  std::unordered_map<const ZEvent *, std::pair<unsigned, unsigned>> event_to_position; ////
 public:
  bool hasEvent(const ZEvent *ev) const;
  bool hasEvent(const ZEventID& ev_id) const;

  // PO
 public:
  ZPartialOrder po;
  const ZPartialOrder& getPo() const { return po; }

  // CACHE
  class Cache {
   public:
    // ML -> thr -> thread-ordered memory-writes
    std::unordered_map
      <SymAddrSize, std::unordered_map
      <unsigned, std::vector<const ZEvent *>>> wm;

    // ML -> thr -> thread-ordered unlocks
    std::unordered_map
      <SymAddrSize, std::unordered_map
      <unsigned, std::vector<const ZEvent *>>> unl;

    // Read -> its local buffer-write
    std::unordered_map
      <const ZEvent *, const ZEvent *> readWB;

    bool empty() const {
      return (wm.empty() && unl.empty() && readWB.empty());
    }
  };
 private:
  Cache cache;
 public:
  const Cache& getCache() const { return cache; }

  bool empty() const { return (lines.empty() && po.empty() && cache.empty()); }
  size_t size() const { return lines.size(); }

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  ~ZGraph() {};
  // Empty
  ZGraph() = delete;
  ZGraph(MemoryModel model);
  // Moving
  ZGraph(ZGraph&& oth);

  ZGraph(const ZGraph& oth) = delete;
  ZGraph& operator=(ZGraph&& oth) = delete;
  ZGraph& operator=(const ZGraph& oth) = delete;

  /* *************************** */
  /* GRAPH CONSTRUCTION          */
  /* *************************** */

  bool constructed = false;

  void construct(const std::vector<std::unique_ptr<ZEvent>>& trace,
                 int pre_tau_limit, std::set<int> causes_after);

  void add_reads_from_edges(const ZAnnotation& annotation);

  /* *************************** */
  /* MAIN ALGORITHM              */
  /* *************************** */

  std::set<ZEventID> get_mutations(
    const ZEvent * const ev, const ZEventID base_obs,
    int mutations_only_from_idx) const;

  std::pair<std::set<int>, std::set<const ZEvent *>> get_causes_after(
    const ZEvent * const readlock, const ZEventID& mutation,
    const std::vector<std::unique_ptr<ZEvent>>& tau) const;

  // In thread thr (and specific aux), starting from ev and going back (ev,ev-1,...,1,0),
  // return first memory-write to ml (resp. return its index in cache.wm)
  int getTailWindex(const SymAddrSize& ml, unsigned thr, int evX) const;
  const ZEvent * const getTailW(const SymAddrSize& ml, unsigned thr, int ev) const;

  // In thread thr (and specific aux), from all memory-writes
  // conflicting with read that do not happen after the read in partial,
  // return the latest one (resp. return its index in cache.wm)
  int getLatestNotAfterIndex(const ZEvent *read, unsigned thr, const ZPartialOrder& partial) const;
  const ZEvent * const getLatestNotAfter(const ZEvent *read, unsigned thr, const ZPartialOrder& partial) const;

  // Get the local buffer-write of the read
  // Returns nullptr if there is no local buffer-write before read
  const ZEvent * const getLocalBufferW(const ZEvent *read) const;

  std::string to_string() const { return po.to_string(); }
  void dump() const { po.dump(); }
};
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZGraph& gr);

#endif // __Z_GRAPH_H__
