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
  ZClosure(const ZAnnotation& annotation, ZPartialOrder& partialOrder);
  ZClosure(const ZClosure& oth) = delete;
  ZClosure& operator=(ZClosure& oth) = delete;
  ZClosure(ZClosure&& oth) = delete;
  ZClosure& operator=(ZClosure&& oth) = delete;

 private:
  std::pair<const ZEvent *, const ZEvent *>
    getObs(const ZEventID& id);
  bool ruleOne
    (const ZEvent *ev, const ZEventID& obs);
  std::pair<bool, bool> ruleTwo
    (const ZEvent *read, const ZEventID& obs);
  std::pair<bool, bool> ruleThree
    (const ZEvent *read, const ZEventID& obs);
  std::pair<bool, bool> rulesTwoThree
    (const ZEvent *read, const ZEventID& obs);

 public:
  bool close();

  const ZAnnotation& an;
  const ZGraph& gr;
  ZPartialOrder& po;

  unsigned added_edges = 0;
  unsigned iterations = 0;
};

#endif // __Z_CLOSURE_H__
