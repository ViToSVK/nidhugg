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

#ifndef _VC_ANNOTATIONNEG_H_
#define _VC_ANNOTATIONNEG_H_

#include "VCAnnotation.h"

#include <vector>

class VCAnnotationNeg {
 public:
  VCAnnotationNeg()
    : forbids_init(false) {}

  VCAnnotationNeg(const VCAnnotationNeg& oth) = default;
  VCAnnotationNeg(VCAnnotationNeg&& oth) = default;

  bool forbidsInitialEvent(const Node * readnd) const {
    assert(isRead(readnd->getEvent()));
    auto key = VCIID(readnd->getProcessID(), readnd->getEventID());
    return forbids_init.count(key);
  }

  bool forbids(const Node * readnd, const Node * writend) const {
    assert(isRead(readnd->getEvent()));
    assert(writend->getEvent() &&
           "should call the special function for init event");
    auto key = VCIID(readnd->getProcessID(), readnd->getEventID());
    auto it = mapping.find(key);
    if (it == mapping.end())
      return false;
    if (writend->getProcessID() >= it->second.size())
      return false;

    return (it->second[writend->getProcessID()]
            >= writend->getEventID());
  }

  void update(const Node * readnd, const std::vector<unsigned>& newneg) {
    assert(isRead(readnd->getEvent()));
    auto key = VCIID(readnd->getProcessID(), readnd->getEventID());
    forbids_init.insert(key);
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

 private:
  std::unordered_set<VCIID> forbids_init;
  std::unordered_map<VCIID, std::vector<unsigned>> mapping;

};

#endif // _VC_ANNOTATIONNEG_H_
