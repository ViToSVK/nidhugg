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


class ZAnnotation {
 public:
  using MappingT = std::map<ZEventID, ZEventID>;
  using LastLockT = std::unordered_map<SymAddrSize, ZEventID>;

 private:
  MappingT mapping;
  LastLockT lastlock;

 public:
  ZAnnotation() = default;
  ZAnnotation(const ZAnnotation&) = default;
  ZAnnotation(ZAnnotation&& a) = default;
  ZAnnotation& operator=(const ZAnnotation&) = delete;
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

  void add(const ZEventID& ev_id, const ZEventID& obs_id);
  void add(const ZEvent *ev, const ZEvent *obs);

  bool defines(const ZEventID& ev_id) const;
  bool defines(const ZEvent *ev) const;

  const ZEventID& obs(const ZEventID& ev_id) const;
  const ZEventID& obs(const ZEvent *ev) const;

  /* *************************** */
  /* LAST LOCK                   */
  /* *************************** */

  void set_last_lock(const ZEvent *ev);

  const ZEventID& last_lock(const ZEvent *ev) const;

  bool is_last_lock(const ZEvent *ev) const;

  bool location_has_some_lock(const ZEvent *ev) const;
};

#endif // _Z_ANNOTATION_H_
