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

#ifndef _Z_ANNOTATION_H_
#define _Z_ANNOTATION_H_

#include <unordered_map>
#include <unordered_set>
#include <sstream>

#include "ZEvent.h"
#include "ZHelpers.h"


class ZObs {
 public:
  ZObs() = delete;
  ZObs(unsigned thread_idx, unsigned event_idx)
    : thr(thread_idx), ev(event_idx) {}

  ZObs(const ZObs& oth) = default;
  ZObs(ZObs&& oth) = default;
  ZObs& operator=(const ZObs& oth) = delete;
  ZObs& operator=(ZObs&& oth) = delete;

  const unsigned thr;
  const unsigned ev;

  std::string to_string() const {
    std::stringstream res;
    res << "[" << thr << ",-1," << ev << ",]";
    return res.str();
  }

  bool operator==(const ZObs& oth) const {
    return (thr == oth.thr && ev == oth.ev);
  }
};


namespace std {
  template <>
  struct hash<ZObs> {
    std::size_t operator()(const ZObs& obs) const {
      return (hash<unsigned>()(obs.thr) << 12) +
             hash<unsigned>()(obs.ev);
    }
  };
}


class ZAnnotation {
 public:

  using MappingT = std::unordered_map<ZObs, ZObs>;
  using LastLockT = std::unordered_map<SymAddrSize, ZObs>;

 private:

  MappingT mapping;
  LastLockT lastlock;

 public:

  ZAnnotation() = default;
  ZAnnotation(const ZAnnotation&) = default;
  ZAnnotation& operator=(const ZAnnotation&) = delete;
  ZAnnotation(ZAnnotation&& a) = default;
  ZAnnotation& operator=(ZAnnotation&& a) = delete;

  bool empty() const { return mapping.empty(); }
  size_t size() const { return mapping.size(); }

  typedef MappingT::iterator iterator;
  typedef MappingT::const_iterator const_iterator;
  iterator begin() { return mapping.begin(); }
  const_iterator begin() const { return mapping.begin(); }
  iterator end() { return mapping.end(); }
  const_iterator end() const { return mapping.end(); }

  void dump() const;


  /* *************************** */
  /* MAPPING                     */
  /* *************************** */

  void add(const ZEvent *ev, const ZObs& obs);

  void add(const ZEvent *ev, const ZEvent *obsEv);

  bool defines(unsigned thrid, unsigned evid) const;

  bool defines(const ZEvent *ev) const;

  const ZObs& getObs(unsigned thrid, unsigned evid) const;

  const ZObs& getObs(const ZEvent *ev) const;


  /* *************************** */
  /* LAST LOCK                   */
  /* *************************** */

  void setLastLock(const ZEvent *ev);

  const ZObs& getLastLock(const ZEvent *ev) const;

  bool isLastLock(const ZEvent *ev) const;

  bool locationHasSomeLock(const ZEvent *ev) const;
};

#endif // _Z_ANNOTATION_H_
