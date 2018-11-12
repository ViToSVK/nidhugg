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

#ifndef _VC_VALCLOSURE_H
#define _VC_VALCLOSURE_H

#include "Debug.h"
#include "VCGraphVclock.h"

class VCValClosure {
 public:

  VCValClosure(const VCGraphVclock& gr, const VCAnnotation& an)
  : graph(gr), annotation(an) {}

  VCValClosure(const VCValClosure& oth) = default;
  VCValClosure& operator=(VCValClosure& oth) = delete;
  VCValClosure(VCValClosure&& oth) = default;
  VCValClosure& operator=(VCValClosure&& oth) = delete;

 private:

  inline bool isGood(const Node * writend, const VCAnnotation::Ann& ann) {
    assert(writend && "Checking nullptr for good");
    assert(!writend->getEvent() || isWrite(writend->getEvent()));
    auto key = std::pair<unsigned, unsigned>(writend->getProcessID(),
                                             writend->getEventID());
    if (ann.goodLocal && *(ann.goodLocal) == key) {
      assert(ann.loc != VCAnnotation::Loc::REMOTE);
      return true;
    }
    if (ann.goodRemote.count(key)) {
      assert(ann.loc != VCAnnotation::Loc::LOCAL);
      return true;
    }
    return false;
  }

  void prepare(const PartialOrder& po, const Node * newread);

  void prepareBounds(const PartialOrder& po, const Node * readnd);

  void updateBounds(const PartialOrder& po, const Node * readnd);

  std::pair<bool, bool> ruleOne
    (const PartialOrder& po, const Node * readnd, const VCAnnotation::Ann& ann);

  std::pair<bool, bool> ruleTwo
    (const PartialOrder& po, const Node * readnd, const VCAnnotation::Ann& ann);

  std::pair<bool, bool> ruleThree
    (const PartialOrder& po, const Node * readnd, const VCAnnotation::Ann& ann);

  std::pair<bool, bool> rules
    (const PartialOrder& po, const Node * readnd, const VCAnnotation::Ann& ann);

 public:

  void valClose(const PartialOrder& po, const Node * newread,
                const VCAnnotation::Ann * newval);

  void valCloseLock(const PartialOrder& po,
                    const Node * locknode,
                    const Node * lastunlocknode);

  const VCGraphVclock& graph;

  const VCAnnotation& annotation;

  bool closed;

  // Bounds for NONROOT reads
  // on visibility of ROOT writes
  std::unordered_map
  <const Node *, std::pair<int, int>> wBounds;
};

#endif // _VC_VALCLOSURE_H
