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


using LineT = std::vector<const ZEvent *>;
using LinesT = std::vector<LineT>;


class ZGraph {
  friend class ZPartialOrder;
  // LINES (changed at recursion child with new ZEvent pointers from new trace)
 private:
  const ZEvent _init; ////
  LinesT _lines; ////
 public:
  const ZEvent *initial() const { return &_init; }
  const LineT& operator()(const CPid& cpid) const;
  //
  const ZEvent *event(const CPid& cpid, int event_id) const;
  const ZEvent *event(const ZEventID& id) const;
  const ZEvent *unlock_of_this_lock(const ZEventID& id) const;
  void add_line(const ZEvent *ev);
  void add_event(const ZEvent *ev);
  void replace_event(const ZEvent *oldEv, const ZEvent *newEv);
  bool has_event(const CPid& cpid, int event_id) const;
  bool has_event(const ZEvent *ev) const;
  std::map<CPid, int> thread_sizes_minus_one() const;
  void shrink();

  // CPID->LINE_ID (retained accross recursion children)
 private:
  std::unordered_map<CPid, unsigned> _cpid_to_line; ////
  std::vector<CPid> _line_to_cpid; ////
  unsigned line_id(const CPid& cpid) const;
  unsigned line_id(const ZEvent *ev) const;
 public:
  const CPid& line_id_to_cpid(unsigned line_id) const;
  std::set<CPid> threads() const;
  bool has_thread(const CPid& cpid) const;

  // PO
 private:
  ZPartialOrder _po; ////
 public:
  ZPartialOrder copy_po() const { return ZPartialOrder(_po); }

  // CACHE
  class Cache {
   public:
    // ML -> CPid -> thread-ordered writes
    std::unordered_map
      <SymAddrSize, std::unordered_map
      <CPid, std::vector<const ZEvent *>>> writes;
    // Read -> its local buffer-write
    std::unordered_map
      <const ZEvent *, const ZEvent *> local_write;
    bool empty() const {
      return (writes.empty() && local_write.empty());
    }
  };
 private:
  Cache _cache; ////
 public:
  const Cache& cache() const { return _cache; }

  bool empty() const {
    return (_lines.empty() && _cpid_to_line.empty() && _line_to_cpid.empty() &&
            _po.empty() && _cache.empty());
  }
  size_t size() const { return _lines.size(); }
  size_t events_size() const {
    size_t res = 0;
    for (auto& ln : _lines)
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
  ZGraph(const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& trace);
  // Moving
  ZGraph(ZGraph&& oth) = default;

  // Extending
  // Partial order that will be moved
  // Trace and annotation that will extend this copy of the graph
  ZGraph(const ZGraph& oth,
         ZPartialOrder&& po,
         const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& trace,
         const ZAnnotation& annotation);

  ZGraph(const ZGraph& oth) = delete;
  ZGraph& operator=(ZGraph&& oth) = delete;
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
  void trace_to_po(const std::vector<std::unique_ptr<ZEvent>>& trace, const ZAnnotation *annotation_ptr);


  /* *************************** */
  /* MAIN ALGORITHM              */
  /* *************************** */

 public:

  // In thread cpid, starting from evX and going back (evX,evX-1,...,1,0),
  // return first write to ml (resp. return its index in cache.wm)
  int get_tailw_index(const SymAddrSize& ml, const CPid& cpid, int evX) const;
  const ZEvent *get_tailw(const SymAddrSize& ml, const CPid& cpid, int evX) const;

  // In thread cpid, from all writes conflicting with read
  // that do not happen after the read in po,
  // return the latest one (resp. return its index in cache.wm)
  int get_latest_not_after_index(const ZPartialOrder& po, const ZEvent *read, const CPid& cpid) const;
  const ZEvent *get_latest_not_after(const ZPartialOrder& po, const ZEvent *read, const CPid& cpid) const;

  // Get the local write of the read
  // Returns nullptr if there is no local write before read
  const ZEvent *get_local_write(const ZEvent *read) const;

  // Returns events-to-mutate in a specified order
  std::list<const ZEvent *> events_to_mutate(const ZAnnotation& annotation) const;

  // Collect all write events visible to the read
  std::set<const ZEvent *> mutation_candidates_collect
  (const ZPartialOrder& po, const ZEvent *read,
   const std::set<ZEventID>& check_if_any_is_visible) const;

  // Filter out candidates forbidden to read by negative annotation
  void mutation_candidates_filter_by_negative
  (const ZEvent *read, std::set<const ZEvent *>& candidates,
   const ZAnnotationNeg& negative) const;

  // Returns mutation candidates for a read node grouped by value
  std::set<ZAnn> mutation_candidates_grouped
  (const ZEvent *read, const ZAnnotationNeg& negative) const;

  std::string to_string() const { return _po.to_string(); }
  void dump() const { _po.dump(); }
};
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZGraph& gr);

#endif // __Z_GRAPH_H__
