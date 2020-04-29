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

#ifndef _Z_ANNOTATION_H_
#define _Z_ANNOTATION_H_

#include <unordered_map>
#include <unordered_set>
#include <sstream>

#include "ZEvent.h"
#include "ZHelpers.h"


class ZAnn {
 public:
  ZAnn(int value, const std::set<ZEventID>& goodwrites)
  : value(value), goodwrites(goodwrites) { assert(!goodwrites.empty()); }

  ZAnn() = delete;
  ZAnn(const ZAnn& a) = default;
  ZAnn(ZAnn&& a) = default;
  ZAnn& operator=(const ZAnn& a) = delete;
  ZAnn& operator=(ZAnn&& a) = delete;

  int value;
  std::set<ZEventID> goodwrites;

  std::string to_string() const;
  void dump() const;

  bool operator==(const ZAnn& c) const { return (compare(c) == 0); }
  bool operator!=(const ZAnn& c) const { return (compare(c) != 0); }
  bool operator<(const ZAnn& c) const { return (compare(c) < 0); }
  bool operator<=(const ZAnn& c) const { return (compare(c) <= 0); }
  bool operator>(const ZAnn& c) const { return (compare(c) > 0); }
  bool operator>=(const ZAnn& c) const { return (compare(c) >= 0); }
  int compare(const ZAnn& c) const;
};
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZAnn& ann);

class ZAnnotation {
 public:
  using MappingT = std::map<ZEventID, ZAnn>;
  using LastLockT = std::unordered_map<SymAddrSize, ZEventID>;

 private:
  MappingT mapping;
  LastLockT lastlock;

 public:
  ZAnnotation() = default;
  ZAnnotation(const ZAnnotation& a) = default;
  ZAnnotation(ZAnnotation&& a) = default;
  ZAnnotation& operator=(const ZAnnotation& a) = delete;
  ZAnnotation& operator=(ZAnnotation&& a) = delete;

  bool empty() const { return mapping.empty() && lastlock.empty(); }
  size_t size() const { return mapping.size(); }

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

  void add(const ZEventID& ev_id, const ZAnn& ann);

  bool defines(const ZEventID& ev_id) const;
  bool defines(const ZEvent *ev) const;

  const ZAnn& ann(const ZEventID& ev_id) const;
  const ZAnn& ann(const ZEvent *ev) const;

  /* *************************** */
  /* LAST LOCK                   */
  /* *************************** */

  void set_last_lock(const ZEvent *ev);

  const ZEventID& last_lock(const ZEvent *ev) const;

  bool is_last_lock(const ZEvent *ev) const;

  bool location_has_some_lock(const ZEvent *ev) const;
};
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZAnnotation& annot);

#endif // _Z_ANNOTATION_H_
