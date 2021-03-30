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

#ifndef _Z_PARTIAL_ORDER_H_
#define _Z_PARTIAL_ORDER_H_

#include <vector>

#include "ZAnnotationNeg.h"

class ZExplorer;

class ZGraph;
class ZPartialOrder {
 public:
  // ZGRAPH REFERENCE (given at constructor-time)
  const ZGraph& graph;  ////
  static ZGraph graphDummy;

  // SUCCESSOR PREDECESSOR
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
  using ClockT = std::vector<std::vector<std::vector<int>>>;
 private:
  ClockT _succ;
  ClockT _pred;
  // For each thread, up until what index (not including that index)
  // the events of the corresponding thread are closure-safe
  std::vector<int> _closure_safe_until;
  // How many events of the given line does the partial order have?
  std::vector<int> _line_sizes;
  int line_size(unsigned line_id) const;
  // Iterate through threads that this partial order spans
  std::vector<CPid> _threads_spanned;
 public:
  const std::vector<CPid>& threads_spanned() const;
  std::map<CPid, int> thread_sizes_minus_one() const;
  // Does the partial order span this thread?
  bool spans_thread(const CPid& cpid) const;
  // How many events of the given thread does the partial order have?
  int thread_size(const CPid& cpid) const;
  // Does the partial order span this event?
  bool spans_event(const ZEvent *ev) const;
  // Is this annotated read provably closure-safe?
  bool is_closure_safe(const ZEvent *read) const;
  // Smallest (i.e. earliest) successor. Returns:
  // Event pointer (nullptr if no successor)
  // Event id (INT_MAX if no successor)
  std::pair<const ZEvent *, int> succ(const ZEvent *from, const CPid& to_cpid) const;
  // Biggest (i.e. latest) predecessor. Returns:
  // Event pointer (nullptr if no predecessor)
  // Event id (-1 if no predecessor)
  std::pair<const ZEvent *, int> pred(const ZEvent *to, const CPid& from_cpid) const;
  // Is there an edge from -> to ?
  bool has_edge(const ZEvent *from, const ZEvent *to) const;
  // Are ev1 and ev2 ordered ?
  bool are_ordered(const ZEvent *ev1, const ZEvent *ev2) const;
  // Add a new edge from -> to (and all edges transitively following)
  void add_edge(const ZEvent *from, const ZEvent *to);
 private:
  // Helper function to maintain transitivity in the partial order
  void add_edge_help(const ZEvent *from, const ZEvent *to);
 public:
  // When creating PO from trace
  void add_line(const CPid& cpid);
  void add_line(const ZEvent * ev);
  void inherit_lines(const ZPartialOrder& oth);
  void add_event(const ZEvent * ev);
  void shrink();
  // Extending when reusing trace and graph
  void extend(const ZEvent *read_lock, const ZAnnotation& mutated_annotation,
              const ZExplorer& explorer, const ZPartialOrder& po_full);
  // Extending but nothing more to annotate, only check backtrack points
  void process_remaining_events_for_backtrack_points
  (const ZExplorer& explorer, const ZPartialOrder& po_full) const;
  // Is there something remaining to annotate?
  bool something_to_annotate(const ZAnnotation& annotation) const;

  // Empty
  ZPartialOrder();
  // Initial
  ZPartialOrder(const ZGraph& graph);
  // When extending
  ZPartialOrder(ZPartialOrder&& oth, const ZGraph& graph);

  ZPartialOrder(ZPartialOrder&& oth) = default;
  ZPartialOrder(const ZPartialOrder& oth) = default;
  ZPartialOrder& operator=(ZPartialOrder&& oth) = delete;
  ZPartialOrder& operator=(const ZPartialOrder& oth) = delete;

  // Get po_full restricted only to events of po_part
  ZPartialOrder(const ZPartialOrder& po_part,
                const ZPartialOrder& po_full);

  bool empty() const {
    assert(_succ.size() == _pred.size());
    return _succ.empty() && _closure_safe_until.empty() && _line_sizes.empty();
  }
  size_t size() const {
    assert(_succ.size() == _pred.size());
    return _succ.size();
  }
  size_t events_size() const {
    size_t res = 0;
    for (const CPid& cpid : threads_spanned())
      res += thread_size(cpid);
    return res;
  }

  std::string to_string() const;
  void dump() const;
};
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZPartialOrder& po);


class POcomp {
  const ZPartialOrder& po;
 public:
  POcomp(const ZPartialOrder& po) : po(po) {}
  bool operator() (const ZEvent *from, const ZEvent *to) const {
    return (po.has_edge(from, to));
  }
};

#endif // _Z_PARTIAL_ORDER_H_
