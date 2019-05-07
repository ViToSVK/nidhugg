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

#include "ZClosure.h"

/* *************************** */
/* RULE 1                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleOne
(const PartialOrder& po, const Node * readnd,
 const ZAnnotation::Ann& ann)
{
  assert(isRead(readnd));
  // Rule1: r observes w ...
  // make w -> r
  const Node *good = getGood(ann);
  if (graph.hasEdge(good, readnd, po))
    return {false, false}; // done, no change
  if (graph.hasEdge(readnd, good, po))
    return {true, false}; // impossible
  graph.addEdge(good, readnd, po);
  return {false, true}; // done, change
}

/* *************************** */
/* RULE 2                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleTwo
(const PartialOrder& po, const Node * readnd,
 const ZAnnotation::Ann& ann)
{
  assert(isRead(readnd));
  // Rule2: r observes w ...
  // w -> w' implies r -> w'
  const Node *good = getGood(ann);
  bool change = false;

  while (true) {
    auto tails = graph.getTailWrites(readnd, po);
    const Node *roottail = tails.first;
    const auto& nonroottails = tails.second;

    // Is Good a tail?
    if (roottail && roottail == good)
      break;
    bool found = false;
    for (const Node * nonroottail : nonroottails)
      if (nonroottail && nonroottail == good) {
        found = true;
        break;
      }
    if (found) break;

    // Good is not a tail, for all tails
    // such that w -> w' make r -> w'
    if (roottail && graph.hasEdge(good, roottail, po)) {
      if (graph.hasEdge(roottail, readnd, po))
        return {true, false}; // impossible
      assert(!graph.hasEdge(readnd, roottail, po));
      graph.addEdge(readnd, roottail, po);
      change = true;
    }
    for (const Node * nonroottail : nonroottails)
      if (nonroottail && graph.hasEdge(good, nonroottail, po)) {
        if (graph.hasEdge(nonroottail, readnd, po))
          return {true, false}; // impossible
        assert(!graph.hasEdge(readnd, nonroottail, po));
        graph.addEdge(readnd, nonroottail, po);
        change = true;
      }
  }

  return {false, change}; // done, change-bool
}

/* *************************** */
/* RULE 3                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleThree
(const PartialOrder& po, const Node * readnd,
 const ZAnnotation::Ann& ann)
{
  assert(isRead(readnd));
  // Rule3: r observes w ...
  // w' -> r implies w' -> w
  const Node *good = getGood(ann);
  bool change = false;

  auto heads = graph.getHeadWrites(readnd, po);
  const Node *roothead = heads.first;
  const auto& nonrootheads = heads.second;

  bool found = false;
  std::unordered_set<const Node *> addEdgeToGood;
  if (roothead) {
    if (roothead == good)
      found = true;
    else {
      assert(!graph.hasEdge(readnd, roothead, po));
      if (graph.hasEdge(roothead, readnd, po))
        addEdgeToGood.insert(roothead);
    }
  }
  for (const Node * nonroothead : nonrootheads) {
    if (nonroothead) {
      if (nonroothead == good)
        found = true;
      else {
        assert(!graph.hasEdge(readnd, nonroothead, po));
        if (graph.hasEdge(nonroothead, readnd, po))
          addEdgeToGood.insert(nonroothead);
      }
    }
  }

  if (!found)
    return {true, false}; // impossible

  for (const Node *badhead : addEdgeToGood) {
    assert(!graph.areOrdered(badhead, good, po));
    graph.addEdge(badhead, good, po);
    change = true;
  }

  #ifndef NDEBUG
  heads = graph.getHeadWrites(readnd, po);
  roothead = heads.first;
  found = false;
  if (roothead && roothead == good)
    found = true;
  else
    assert(!roothead || !graph.hasEdge(roothead, readnd, po));
  for (const Node * nonroothead : heads.second) {
    if (nonroothead && nonroothead == good)
      found = true;
    else
      assert(!nonroothead || !graph.hasEdge(nonroothead, readnd, po));
  }
  assert(found);
  #endif

  return {false, change}; // done, change-bool
}

/* *************************** */
/* RULES                       */
/* *************************** */

std::pair<bool, bool> ZClosure::rules
(const PartialOrder& po, const Node * readnd,
 const ZAnnotation::Ann& ann)
{
  assert(readnd && isRead(readnd));

  bool change = false;
  // Rule1 is done only the first time
  // Rule2
  auto res = ruleTwo(po, readnd, ann);
  if (res.first) return {true, false};
  if (res.second) change = true;
  //Rule3
  res = ruleThree(po, readnd, ann);
  if (res.first) return {true, false};
  if (res.second) change = true;

  return {false, change};
}


/* *************************** */
/* VAL CLOSE                   */
/* *************************** */

void ZClosure::valClose
(const PartialOrder& po, const Node * newread,
 const ZAnnotation::Ann * newann)
{
  // Rule1 for new annotation
  if (newread) {
    auto res = ruleOne(po, newread, *newann);
    if (res.first) { closed = false; return; }
  }

  bool change = true;
  while (change) {
    change = false;
    for (const auto& key_ann : annotation) {
      assert(graph.closureSafeUntil.size() > key_ann.first.first);
      if (graph.closureSafeUntil[key_ann.first.first]
          < (int) key_ann.first.second) {
        auto res = rules(po, graph.getNode(key_ann.first), key_ann.second);
        if (res.first) { closed = false; return; }
        if (res.second) change = true;
      }
    }
    if (newread) {
      auto res = rules(po, newread, *newann);
      if (res.first) { closed = false; return; }
      if (res.second) change = true;
    }
  }

  assert(!change);
  closed = true;
}

/* *************************** */
/* VAL CLOSE LOCK              */
/* *************************** */

void ZClosure::valCloseLock(const PartialOrder& po,
                                const Node * locknode,
                                const Node * lastunlocknode)
{
  assert(isLock(locknode) && isUnlock(lastunlocknode) &&
         sameMl(locknode, lastunlocknode));

  if (graph.hasEdge(locknode, lastunlocknode, po)) {
    closed = false;
    return;
  } else if (graph.hasEdge(lastunlocknode, locknode, po)) {
    closed = true;
    return;
  }

  graph.addEdge(lastunlocknode, locknode, po);

  valClose(po, nullptr, nullptr);
}
