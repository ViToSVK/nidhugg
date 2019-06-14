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

#ifndef _Z_ANNOTATIONNEG_H_
#define _Z_ANNOTATIONNEG_H_

#include <vector>

#include "ZAnnotation.h"


class ZAnnotationNeg {
 public:
  ZAnnotationNeg() = default;
  ZAnnotationNeg(const ZAnnotationNeg& oth) = default;
  ZAnnotationNeg(ZAnnotationNeg&& oth) = default;
  ZAnnotationNeg& operator=(const ZAnnotationNeg&) = delete;
  ZAnnotationNeg& operator=(ZAnnotationNeg&& a) = delete;

  bool forbidsInitialEvent(const ZEvent *readev) const;

  bool forbids(const ZEvent *readev, const ZEvent *writeev) const;

  void update(const ZEvent *readev, std::vector<unsigned>&& newneg);

  bool empty() const { return mapping.empty(); }

  std::string to_string() const;
  void dump() const;

  using MappingT = std::map<ZObs, std::vector<unsigned>>;
  typedef MappingT::iterator iterator;
  typedef MappingT::const_iterator const_iterator;

  iterator begin() { return mapping.begin(); }
  const_iterator begin() const { return mapping.begin(); }
  iterator end() { return mapping.end(); }
  const_iterator end() const { return mapping.end(); }

 private:
  MappingT mapping;

};

#endif // _Z_ANNOTATIONNEG_H_
