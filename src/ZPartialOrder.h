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

#ifndef _Z_PARTIAL_ORDER_H_
#define _Z_PARTIAL_ORDER_H_

#include "ZBasis.h"


class ZPartialOrder {
 public:
  // ZBASIS REFERENCE (given at constructor-time)
  const ZBasis& basis;  ////
  typedef std::vector<std::vector<std::vector<int>>> ClockT;


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
 private:
  ClockT _succ;
  ClockT _pred;
  std::pair<const ZEvent *, int> succ(const ZEvent *from, unsigned to_line) const;
  std::pair<const ZEvent *, int> pred(const ZEvent *to, unsigned to_line) const;
 public:
  std::pair<const ZEvent *, int> succ(const ZEvent *from, unsigned to_thread, int to_aux) const;
  std::pair<const ZEvent *, int> pred(const ZEvent *to, unsigned to_thread, int to_aux) const;
  bool hasEdge(const ZEvent *from, const ZEvent *to) const;
  bool areOrdered(const ZEvent *ev1, const ZEvent *ev2) const;
  void addEdge(const ZEvent *from, const ZEvent *to);
 private:
  void addEdgeHelp(const ZEvent *from, const ZEvent *to);
 public:
  // When creating PO from trace
  void addLine(const ZEvent * ev);
  void addEvent(const ZEvent * ev);


 public:
  // ZPartialOrder is not responsible for any resources
  ZPartialOrder() = delete;
  ZPartialOrder(const ZBasis& basis)
    : basis(basis),
    _succ(),
    _pred()
      {}

  ZPartialOrder(const ZPartialOrder& oth, const ZBasis& basis)
    : basis(basis),
    _succ(oth._succ),
    _pred(oth._pred)
      {}

  ZPartialOrder(ZPartialOrder&& oth) = default;
  ZPartialOrder& operator=(ZPartialOrder&& oth) = delete;
  ZPartialOrder(const ZPartialOrder& oth) = delete;
  ZPartialOrder& operator=(const ZPartialOrder& oth) = delete;

  bool empty() const {
    assert(_succ.size() == _pred.size());
    return _succ.empty();
  }
  size_t size() const {
    assert(_succ.size() == _pred.size());
    return _succ.size();
  }

  void dump() const;
};


class POcomp {
  const ZPartialOrder& po;
 public:
  POcomp(const ZPartialOrder& po) : po(po) {}
  bool operator() (const ZEvent *from, const ZEvent *to) const {
    return (po.hasEdge(from, to));
  }
};

#endif // _Z_PARTIAL_ORDER_H_
