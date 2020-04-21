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

#ifndef _Z_ANNOTATIONNEG_H_
#define _Z_ANNOTATIONNEG_H_

#include "ZAnnotation.h"


class ZAnnotationNeg {
 public:
  using NegativeT = std::map<CPid, int>;
  using MappingT = std::map<ZEventID, NegativeT>;

 private:
  MappingT mapping;

 public:
  ZAnnotationNeg() = default;
  ZAnnotationNeg(const ZAnnotationNeg& oth) = default;
  ZAnnotationNeg(ZAnnotationNeg&& oth) = default;
  ZAnnotationNeg& operator=(const ZAnnotationNeg&) = delete;
  ZAnnotationNeg& operator=(ZAnnotationNeg&& a) = delete;

  bool empty() const { return mapping.empty(); }

  using iterator = MappingT::iterator;
  using const_iterator = MappingT::const_iterator;
  iterator begin() { return mapping.begin(); }
  const_iterator begin() const { return mapping.begin(); }
  iterator end() { return mapping.end(); }
  const_iterator end() const { return mapping.end(); }

  std::string to_string() const;
  void dump() const;

  /* *************************** */
  /* MAPPING                     */
  /* *************************** */

  bool forbids_initial(const ZEvent *ev) const;

  bool forbids(const ZEvent *ev, const ZEvent *obs) const;

  void update(const ZEvent *ev, NegativeT&& upd);
};

#endif // _Z_ANNOTATIONNEG_H_
