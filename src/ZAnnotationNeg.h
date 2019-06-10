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

  bool forbidsInitialEvent(const ZEvent *readev) const {
    assert(isRead(readev) || isLock(readev));
    auto key = ZObs(readev->threadID(), readev->eventID());
    return mapping.count(key);
  }

  bool forbids(const ZEvent *readev, const ZEvent *writeev) const {
    assert(writeev && !isInitial(writeev) &&
           "Call the special function for init event");
    assert((isRead(readev) && isWriteB(writeev)) ||
           (isLock(readev) && isUnlock(writeev)));
    assert(readev->auxID() == -1 && writeev->auxID() == -1);
    auto key = ZObs(readev->threadID(), readev->eventID());
    auto it = mapping.find(key);
    if (it == mapping.end())
      return false;
    if (writeev->threadID() >= it->second.size())
      return false;

    return (it->second[writeev->threadID()]
            >= writeev->eventID());
  }

  void update(const ZEvent *readev, const std::vector<unsigned>& newneg) {
    assert(isRead(readev) || isLock(readev));
    auto key = ZObs(readev->threadID(), readev->eventID());
    auto it = mapping.find(key);
    if (it == mapping.end())
      mapping.emplace_hint(it, key, newneg);
    else {
      assert(it->second.size() <= newneg.size());
      it->second.reserve(newneg.size());
      for (unsigned i = 0; i < newneg.size(); ++i)
        if (i < it->second.size()) {
          assert(it->second[i] <= newneg[i]);
          it->second[i] = newneg[i];
        } else
          it->second.push_back(newneg[i]);
    }
  }

  bool empty() const { return mapping.empty(); }

 private:
  std::unordered_map<ZObs, std::vector<unsigned>> mapping;

};

#endif // _Z_ANNOTATIONNEG_H_
