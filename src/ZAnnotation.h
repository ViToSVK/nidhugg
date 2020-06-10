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
  // Captures the buffer-write part of the observation
  // Special case: (CPid(),-1) for observing initial event
  using MappingT = std::map<ZEventID, ZEventID>;

 private:
  MappingT read_mapping;
  MappingT lock_mapping;

 public:
  ZAnnotation() = default;
  ZAnnotation(const ZAnnotation&) = default;
  ZAnnotation(ZAnnotation&& a) = default;
  ZAnnotation& operator=(const ZAnnotation&) = delete;
  ZAnnotation& operator=(ZAnnotation&& a) = delete;

  bool empty() const { return read_mapping.empty() && lock_mapping.empty(); }
  size_t size() const { return read_mapping.size(); }
  size_t read_size() const { return read_mapping.size(); }
  size_t lock_size() const { return lock_mapping.size(); }

  using iterator = MappingT::iterator;
  using const_iterator = MappingT::const_iterator;
  iterator begin() { return read_mapping.begin(); }
  const_iterator begin() const { return read_mapping.begin(); }
  iterator end() { return read_mapping.end(); }
  const_iterator end() const { return read_mapping.end(); }
  const_iterator read_begin() const { return read_mapping.begin(); }
  const_iterator read_end() const { return read_mapping.end(); }
  const_iterator lock_begin() const { return lock_mapping.begin(); }
  const_iterator lock_end() const { return lock_mapping.end(); }

  std::string to_string() const;
  void dump() const;

  /* Comparison implements a total order over ZAnnotations. */
  bool operator==(const ZAnnotation &c) const { return compare(c) == 0; };
  bool operator!=(const ZAnnotation &c) const { return compare(c) != 0; };
  bool operator<(const ZAnnotation &c) const { return compare(c) < 0; };
  bool operator<=(const ZAnnotation &c) const { return compare(c) <= 0; };
  bool operator>(const ZAnnotation &c) const { return compare(c) > 0; };
  bool operator>=(const ZAnnotation &c) const { return compare(c) >= 0; };
  int compare(const ZAnnotation &c) const;

  /* *************************** */
  /* READS                       */
  /* *************************** */

  void add(const ZEventID& ev_id, const ZEventID& obs_id);
  void add(const ZEvent *ev, const ZEvent *obs);

  bool defines(const ZEventID& ev_id) const;
  bool defines(const ZEvent *ev) const;

  const ZEventID& obs(const ZEventID& ev_id) const;
  const ZEventID& obs(const ZEvent *ev) const;

  /* *************************** */
  /* LOCKS                       */
  /* *************************** */

  void lock_add(const ZEventID& ev_id, const ZEventID& obs_id);
  void lock_add(const ZEvent *ev, const ZEvent *obs);

  bool lock_defines(const ZEventID& ev_id) const;
  bool lock_defines(const ZEvent *ev) const;

  const ZEventID& lock_obs(const ZEventID& ev_id) const;
  const ZEventID& lock_obs(const ZEvent *ev) const;
};
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZAnnotation& annot);

#endif // _Z_ANNOTATION_H_
