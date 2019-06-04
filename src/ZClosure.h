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

#ifndef __Z_CLOSURE_H__
#define __Z_CLOSURE_H__

#include "ZGraph.h"


class ZClosure {
 public:

  ZClosure(const ZGraph& graph,
           const ZAnnotation& annotation,
           ZPartialOrder& partialOrder)
  : gr(graph), an(annotation), po(partialOrder) {
    assert(&(po.basis.graph) == &graph);
  }

  ZClosure(const ZClosure& oth) = delete;
  ZClosure& operator=(ZClosure& oth) = delete;
  ZClosure(ZClosure&& oth) = delete;
  ZClosure& operator=(ZClosure&& oth) = delete;

 private:

  std::pair<const ZEvent *, const ZEvent *>
    getObs(const ZObs& obs);

  std::pair<bool, bool> ruleOne
    (const ZEvent *read, const ZObs& obs);

  std::pair<bool, bool> ruleTwo
    (const ZEvent *read, const ZObs& obs);

  std::pair<bool, bool> ruleThree
    (const ZEvent *read, const ZObs& obs);

  std::pair<bool, bool> rules
    (const ZEvent *read, const ZObs& obs);

 public:

  bool close
    (const ZEvent *newread);

  void preClose
    (const ZEvent *ev, const ZEvent *obsEv);
  void preClose
    (const ZEvent *ev, const ZObs& obs);

  const ZGraph& gr;

  const ZAnnotation& an;

  ZPartialOrder& po;
};

#endif // __Z_CLOSURE_H__
