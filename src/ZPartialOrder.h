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

#include "ZAnnotation.h"


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
  typedef std::vector<std::vector<std::vector<int>>> ClockT;
 private:
  ClockT _succ;
  ClockT _pred;
  std::pair<const ZEvent *, int> succ(const ZEvent *from, unsigned to_line) const;
  std::pair<const ZEvent *, int> pred(const ZEvent *to, unsigned to_line) const;
 public:
  // Smallest (i.e. earliest) successor. Returns:
  // Event pointer (nullptr if no successor)
  // Event id (INT_MAX if no successor)
  std::pair<const ZEvent *, int> succ(const ZEvent *from, unsigned to_thread, int to_aux) const;
  // Biggest (i.e. latest) predecessor. Returns:
  // Event pointer (nullptr if no predecessor)
  // Event id (-1 if no predecessor)
  std::pair<const ZEvent *, int> pred(const ZEvent *to, unsigned to_thread, int to_aux) const;
  // Is there an edge from -> to ?
  bool hasEdge(const ZEvent *from, const ZEvent *to) const;
  // Are ev1 and ev2 ordered ?
  bool areOrdered(const ZEvent *ev1, const ZEvent *ev2) const;
  // Add a new edge from -> to (and all edges transitively following)
  void addEdge(const ZEvent *from, const ZEvent *to);
 private:
  // Helper function to maintain transitivity in the partial order
  void addEdgeHelp(const ZEvent *from, const ZEvent *to);
 public:
  // When creating PO from trace
  void addLine(const ZEvent * ev);
  void addEvent(const ZEvent * ev);
  void shrink();

 public:
  // Initial
  ZPartialOrder(const ZGraph& graph);
  // Empty
  ZPartialOrder();

  ZPartialOrder(const ZPartialOrder& oth) = default;
  ZPartialOrder(ZPartialOrder&& oth) = default;
  ZPartialOrder& operator=(ZPartialOrder&& oth) = delete;
  ZPartialOrder& operator=(const ZPartialOrder& oth) = delete;

  bool empty() const {
    assert(_succ.size() == _pred.size());
    return _succ.empty();
  }
  size_t size() const {
    assert(_succ.size() == _pred.size());
    return _succ.size();
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
    return (po.hasEdge(from, to));
  }
};

#endif // _Z_PARTIAL_ORDER_H_
