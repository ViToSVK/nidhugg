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

#ifndef __Z_LINEARIZATION_H__
#define __Z_LINEARIZATION_H__

#include "ZGraph.h"


class ZLinearization {
 public:

  ZLinearization(const ZAnnotation& annotation,
                 ZPartialOrder& partialOrder)
  : an(annotation), gr(partialOrder.basis.graph),
    ba(partialOrder.basis), po(partialOrder) {}

  ZLinearization(const ZLinearization& oth) = delete;
  ZLinearization& operator=(ZLinearization& oth) = delete;
  ZLinearization(ZLinearization&& oth) = delete;
  ZLinearization& operator=(ZLinearization&& oth) = delete;

  const ZAnnotation& an;
  const ZGraph& gr;
  const ZBasis& ba;
  const ZPartialOrder& po;

  std::vector<ZEvent> linearizeTSO() const;
  std::vector<ZEvent> linearizePSO() const;
};

#endif // __Z_LINEARIZATION_H__
