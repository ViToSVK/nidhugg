/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2018 Viktor Toman
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

#include "VCBasis.h"
#include "VCHelpers.h"
#include "Debug.h"

#include <unordered_map>
#include <unordered_set>

typedef std::pair<unsigned, unsigned> VCIID;

namespace std {
  template <>
  struct hash<VCIID> {
    std::size_t operator()(const VCIID& vciid) const {
      return (hash<unsigned>()(vciid.first) << 12) +
             hash<unsigned>()(vciid.second);
    }
  };
}

class ZAnnotation {
 public:
  enum class Loc {
    LOCAL, REMOTE, ANY
  };

  class Ann {
   public:
    ~Ann() { delete goodLocal; }

    Ann(int v, Loc l, std::unordered_set<VCIID>&& gr, bool haslocal, VCIID gl)
      : value(v), loc(l), goodRemote(std::move(gr)),
      goodLocal(haslocal?(new VCIID(gl.first, gl.second)):nullptr)
        {}

    Ann(const Ann& oth)
      : value(oth.value),
      loc(oth.loc),
      goodRemote(oth.goodRemote),
      goodLocal(oth.goodLocal ?
                new VCIID(oth.goodLocal->first,
                          oth.goodLocal->second) :
                nullptr)
        {}

    Ann(Ann&& oth)
      : value(oth.value),
      loc(oth.loc),
      goodRemote(std::move(oth.goodRemote)),
      goodLocal(oth.goodLocal)
        { oth.goodLocal = nullptr; }

    const int value;
    const Loc loc;
    const std::unordered_set<VCIID> goodRemote;
    const VCIID * goodLocal;

    void dump() const;
  };

  using MappingT = std::unordered_map<VCIID, Ann>;
  using LastLockT = std::unordered_map<SymAddrSize, VCIID>;
  using EverGoodT = std::unordered_map<SymAddrSize, std::unordered_set<VCIID>>;

 private:

  MappingT mapping;
  LastLockT lastlock;
  EverGoodT everGood;

 public:

  ZAnnotation() = default;
  ZAnnotation(const ZAnnotation&) = default;
  ZAnnotation& operator=(const ZAnnotation&) = default;
  ZAnnotation(ZAnnotation&& a) = default;
  ZAnnotation& operator=(ZAnnotation&& a) = default;

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

 public:

  // Retuns VCIIDs of newly everGood writes in ann
  std::unordered_set<VCIID> add(const Node * nd, const Ann& ann) {
    assert(isRead(nd));
    auto key = VCIID(nd->getProcessID(), nd->getEventID());
    auto it = mapping.find(key);
    assert(it == mapping.end());
    assert((ann.goodRemote.size() == 1 && !ann.goodLocal) ||
           (ann.goodRemote.empty() && ann.goodLocal)); // dc
    mapping.emplace_hint(it, key, ann);

    // Maintain the set of everGood writes
    // Collect and return newly everGood writes
    auto result = std::unordered_set<VCIID>();
    auto mlit = everGood.find(nd->getEvent()->ml);
    if (mlit == everGood.end())
      mlit = everGood.emplace_hint(mlit, nd->getEvent()->ml,
                                   std::unordered_set<VCIID>());
    if (ann.goodLocal) {
      auto neweverg = mlit->second.emplace(ann.goodLocal->first,
                                           ann.goodLocal->second);
      if (neweverg.second)
        result.emplace(ann.goodLocal->first, ann.goodLocal->second);
    }
    for (auto& vciid : ann.goodRemote) {
      auto neweverg = mlit->second.emplace(vciid.first, vciid.second);
      if (neweverg.second)
        result.emplace(vciid.first, vciid.second);
    }
    return result;
  }

  bool defines(const Node * nd) const {
    assert(isRead(nd));
    auto key = VCIID(nd->getProcessID(), nd->getEventID());
    return (mapping.find(key) != mapping.end());
  }

  bool defines(unsigned pid, unsigned evid) const {
    auto key = VCIID(pid, evid);
    return (mapping.find(key) != mapping.end());
  }

  bool isEverGood(const Node * nd) const {
    assert(isWrite(nd));
    auto mlit = everGood.find(nd->getEvent()->ml);
    if (mlit == everGood.end()) {
      return false;
    }
    return mlit->second.count(VCIID(nd->getProcessID(),
                                    nd->getEventID()));
  }

  const Ann& getAnn(const Node * nd) const {
    assert(isRead(nd));
    auto key = VCIID(nd->getProcessID(), nd->getEventID());
    auto it = mapping.find(key);
    assert(it != mapping.end());
    return it->second;
  }

  const Ann& getAnn(unsigned pid, unsigned evid) const {
    auto key = VCIID(pid, evid);
    auto it = mapping.find(key);
    assert(it != mapping.end());
    return it->second;
  }

  /* *************************** */
  /* LAST LOCK                   */
  /* *************************** */

  void setLastLock(const Node * nd) {
    assert(isLock(nd));
    auto it = lastlock.find(nd->getEvent()->ml);
    if (it == lastlock.end())
      lastlock.emplace_hint(it, nd->getEvent()->ml,
                            VCIID(nd->getProcessID(), nd->getEventID()));
    else
      it->second = VCIID(nd->getProcessID(), nd->getEventID());
  }

  std::pair<bool, VCIID> getLastLock(const Node * nd) const {
    assert(isLock(nd));
    auto it = lastlock.find(nd->getEvent()->ml);
    if (it == lastlock.end())
      return {false, VCIID(1337,47)};
    else
      return {true, it->second};
  }

  bool isLastLock(const Node * nd) const {
    assert(isLock(nd));
    auto it = lastlock.find(nd->getEvent()->ml);
    if (it == lastlock.end())
      return false;
    else
      return (it->second == VCIID(nd->getProcessID(), nd->getEventID()));
  }

  bool locationHasSomeLock(const Node * nd) const {
    assert(isLock(nd));
    return (lastlock.find(nd->getEvent()->ml) != lastlock.end());
  }
};

namespace std {
  template <>
  struct hash<std::pair<int, ZAnnotation::Loc>> {
    std::size_t operator()(const std::pair<int, ZAnnotation::Loc>& val_loc) const {
      return hash<int>()(val_loc.first) +
             (size_t) val_loc.second;
    }
  };
}

#endif // _Z_ANNOTATION_H_