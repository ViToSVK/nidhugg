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

#ifndef __Z_CLOSURE_H__
#define __Z_CLOSURE_H__

#include "ZGraph.h"


class ZClosure {
 public:
  ZClosure(const ZAnnotation& annotation,
           ZPartialOrder& partial_order)
  : an(annotation), gr(partial_order.graph), po(partial_order) {}

  ZClosure(const ZClosure& oth) = delete;
  ZClosure(ZClosure&& oth) = delete;
  ZClosure& operator=(ZClosure& oth) = delete;
  ZClosure& operator=(ZClosure&& oth) = delete;

 private:

  std::pair<bool, bool> rule_two
    (const ZEvent *read, const ZEvent *write);

  std::pair<bool, bool> rule_three
    (const ZEvent *read, const ZEvent *write);

  std::pair<bool, bool> rules
    (const ZEvent *read, const ZAnn& ann);

 public:

  bool close_finish
    (const std::map<ZEventID, ZAnn>& reads_with_multiple_good_writes) const;

  bool close
    (const ZEvent *newread);

  void rule_one
    (const ZEvent *ev, const ZAnn& ann);

  void rule_one_lock
    (const ZEvent *lock, const ZEvent *unlock);

  const ZAnnotation& an;
  const ZGraph& gr;
  ZPartialOrder& po;

  unsigned added_edges = 0;
  unsigned iterations = 0;
};

#endif // __Z_CLOSURE_H__
