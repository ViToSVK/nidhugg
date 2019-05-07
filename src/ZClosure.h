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

#ifndef __Z_CLOSURE_H__
#define __Z_CLOSURE_H__

#include "Debug.h"
#include "ZGraph.h"

class ZClosure {
 public:

  ZClosure(const ZGraph& gr, const ZAnnotation& an)
  : graph(gr), annotation(an) {}

  ZClosure(const ZClosure& oth) = default;
  ZClosure& operator=(ZClosure& oth) = delete;
  ZClosure(ZClosure&& oth) = default;
  ZClosure& operator=(ZClosure&& oth) = delete;

 private:

  const Node *getGood(const ZAnnotation::Ann& ann) {
    if (ann.goodLocal) {
      assert(ann.loc != ZAnnotation::Loc::REMOTE);
      if (ann.goodLocal->first == INT_MAX)
        return graph.initial_node;
      auto result = graph.getNode(*(ann.goodLocal));
      assert(isWrite(result));
      return result;
    }
    assert(ann.goodRemote.size() == 1);
    auto result = graph.getNode(*(ann.goodRemote.begin()));
    assert(isWrite(result));
    return result;
  }

  std::pair<bool, bool> ruleOne
    (const PartialOrder& po, const Node * readnd, const ZAnnotation::Ann& ann);

  std::pair<bool, bool> ruleTwo
    (const PartialOrder& po, const Node * readnd, const ZAnnotation::Ann& ann);

  std::pair<bool, bool> ruleThree
    (const PartialOrder& po, const Node * readnd, const ZAnnotation::Ann& ann);

  std::pair<bool, bool> rules
    (const PartialOrder& po, const Node * readnd, const ZAnnotation::Ann& ann);

 public:

  void valClose(const PartialOrder& po, const Node * newread,
                const ZAnnotation::Ann * newval);

  void valCloseLock(const PartialOrder& po,
                    const Node * locknode,
                    const Node * lastunlocknode);

  const ZGraph& graph;

  const ZAnnotation& annotation;

  bool closed;
};

#endif // __Z_CLOSURE_H__
