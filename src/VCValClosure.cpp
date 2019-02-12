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

#include "VCValClosure.h"

void VCValClosure::prepare
(const PartialOrder& po, const Node * newread, const VCAnnotation::Ann * newann)
{
  if (newread) {
    // performing a new annotation
    // caches for all old nonroot annotations are set
    #ifndef NDEBUG
    for (auto& key_ann : annotation)
      assert(key_ann.second.oneGW ||
             graph.getNode(key_ann.first)->getProcessID() == graph.starRoot() ||
             wBounds.count(graph.getNode(key_ann.first)));
    #endif
    assert(newann);
    if (newread->getProcessID() != graph.starRoot() && !newann->oneGW)
      prepareBounds(po, newread);
  } else {
    // either fully preparing an extension according to the
    // previous annotation, or closing after a lock annotation
    for (auto& key_ann : annotation) {
      if (key_ann.second.oneGW)
        continue;
      if (graph.getNode(key_ann.first)->getProcessID() != graph.starRoot())
        prepareBounds(po, graph.getNode(key_ann.first));
    }
  }
}

void VCValClosure::prepareBounds
(const PartialOrder& po, const Node * readnd)
{
  assert(isRead(readnd));
  assert(readnd->getProcessID() != graph.starRoot());

  wBounds.emplace(readnd, std::pair<int, int>
                  (0, graph.wRoot.at(readnd->getEvent()->ml).size()-1));
  updateBounds(po, readnd);
}

void VCValClosure::updateBounds
(const PartialOrder& po, const Node * readnd)
{
  assert(isRead(readnd));
  assert(readnd->getProcessID() != graph.starRoot());
  assert(wBounds.count(readnd));

  int& low = wBounds[readnd].first;
  int& high = wBounds[readnd].second;
  const std::vector<const Node *>&
    wRemote = graph.wRoot.at(readnd->getEvent()->ml);

  // Update wBounds
  // Note: replace with binary searches?
  assert(low >= high || graph.hasEdge(wRemote[low],
                                      wRemote[high], po));
  while (high >= 0) {
    if (graph.hasEdge(readnd, wRemote[high], po))
      --high; // this one is covered (since the read HB it)
    else
      break;
  }
  while (low < (int) wRemote.size()) {
    if (low + 1 < (int) wRemote.size() &&
        graph.hasEdge(wRemote[low + 1], readnd, po))
      ++low; // this one is covered
    else
      break;
  }
  // Assertions
  assert(low <= high+1);
  assert(low >= high || graph.hasEdge(wRemote[low],
                                      wRemote[high], po));
  assert((high == -1 || !graph.hasEdge(readnd, wRemote[high], po))
         && high + 1 <= (int) wRemote.size() &&
         (high + 1 == (int) wRemote.size() ||
          graph.hasEdge(readnd, wRemote[high + 1], po))
         && "badly computed high index into root writes");
  // Our 'low' might be not entirely correct here, as the write on
  // 'low' might be covered from readnd by some nonroot write, but this
  // can't be cheaply checked during update, so we check during actual rules
  assert((low == (int) wRemote.size()
          ||
         (
          graph.hasEdge(wRemote[low], readnd, po)
          &&
          (low + 1 == (int) wRemote.size() ||
           !graph.hasEdge(wRemote[low + 1], readnd, po))
         )
          ||
         (
          !graph.hasEdge(wRemote[low], readnd, po)
          && low == 0
         )) && "badly computed low index into root writes");
}

/* *************************** */
/* RULE 1                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> VCValClosure::ruleOne
(const PartialOrder& po, const Node * readnd,
 const VCAnnotation::Ann& ann)
{
  assert(isRead(readnd));
  assert(readnd->getProcessID() == graph.starRoot() || ann.oneGW ||
         wBounds.count(readnd));
  // Rule1: there exists a write such that
  // (i) GOOD (ii) HEAD (iii) HB read
  //
  // Note: a local write is always
  // either *HEAD* or *covered*

  if (ann.oneGW && ann.loc != VCAnnotation::Loc::LOCAL) {
    /* ***************** */
    /* DATACENTRIC       */
    /* ***************** */
    // Rule1: r observes w ...
    // make w -> r
    assert(readnd->getProcessID() == graph.starRoot() &&
           "Rule1 done only for root read nodes");
    const Node *good = getGood(ann);

    bool change = false;
    if (graph.hasEdge(readnd, good, po))
      return {true, false}; // impossible
    if (!graph.hasEdge(good, readnd, po)) {
      graph.addEdge(good, readnd, po);
      change = true;
    }

    auto heads = graph.getHeadWrites(readnd, po);
    const Node *roothead = heads.first;
    if (roothead && roothead == good)
      return {false, change}; // done, change bool
    if(heads.second.size() == 1 &&
       *(heads.second.begin()) == good)
      return {false, change}; // done, change bool
    return {true, false}; // impossible
  }

  assert(!ann.oneGW || ann.loc == VCAnnotation::Loc::LOCAL);
  assert(readnd->getProcessID() == graph.starRoot() &&
         "Leaf doesn't need Rule1");

  if (ann.loc == VCAnnotation::Loc::LOCAL) {
    /* ***************** */
    /* ROOT   LOCAL      */
    /* ***************** */

    // Check root head, that's the only way to satisfy
    auto heads = graph.getHeadWrites(readnd, po);
    const Node *roothead = heads.first;
    if (!roothead || !isGood(roothead, ann))
      return {true, false}; // impossible

    assert(graph.hasEdge(roothead, readnd, po));
    return {false, false}; // done, no change
  }

  assert(ann.loc == VCAnnotation::Loc::REMOTE);

  /* ***************** */
  /* ROOT   REMOTE     */
  /* ***************** */

  // Check nonroot head(s)
  auto heads = graph.getHeadWrites(readnd, po);
  const auto& nonrootheads = heads.second;
  if (nonrootheads.size() == 1 &&
      isGood(*(nonrootheads.begin()), ann)) {
    const Node *nonrhead = *(nonrootheads.begin());
    // We have good remote head, now we
    // just need to check the HB edge
    assert(!graph.hasEdge(readnd, nonrhead, po));
    if (graph.hasEdge(nonrhead, readnd, po)) {
      // Edge is already there, we are done
      return {false, false}; // done, no change
    } else {
      // We have to add the edge
      graph.addEdge(nonrhead, readnd, po);
      // This edge doesn't change heads
      #ifndef NDEBUG
      auto newheads = graph.getHeadWrites(readnd, po);
      assert(newheads.first == heads.first);
      assert(newheads.second.size() == 1 &&
             *(newheads.second.begin()) == nonrhead);
      #endif
      return {false, true}; // done, change
    }
  }

  // Either no nonroot head(s) or bad nonroot head(s)
  #ifndef NDEBUG
  if (nonrootheads.size() > 0)
    for (const Node *badNonrootHead : nonrootheads) {
      assert(!isGood(badNonrootHead, ann));
    }
  #endif
  // Remote observable good writes are fully ordered, quicksort them
  auto remoteGoodWrites = std::vector<const Node *>();
  auto itunord = graph.wNonrootUnord.find(readnd->getEvent()->ml);
  assert(itunord != graph.wNonrootUnord.end());
  for (auto& remoteWrite : itunord->second)
    if (isGood(remoteWrite, ann) && graph.isObservableBy(remoteWrite, readnd, po))
      remoteGoodWrites.push_back(remoteWrite);

  if (remoteGoodWrites.size() == 0)
    return {true, false}; // impossible
  #ifndef NDEBUG
  for (auto it1 = remoteGoodWrites.begin();
       it1 != remoteGoodWrites.end(); ++it1)
    for (auto it2 = remoteGoodWrites.begin();
         it2 != it1; ++it2) {
      assert(graph.areOrdered(*it1, *it2, po));
    }
  #endif
  std::sort(remoteGoodWrites.begin(), remoteGoodWrites.end(),
            VCGraphVclock::POcomp(graph, po));

  unsigned rgw_id = 0;
  while (rgw_id < remoteGoodWrites.size()) {
    // We are trying with the write in position rgw_id
    const Node * rgw_current = remoteGoodWrites[rgw_id];
    if (graph.hasEdge(readnd, rgw_current, po)) {
      // Previous all didn't work and this one already HA readnd
      return {true, false}; // impossible
    }
    if (graph.hasEdge(rgw_current, readnd, po)) {
      // This one is covered, we have to go further
      #ifndef NDEBUG
      if (nonrootheads.size() == 0)
        assert(heads.first && graph.hasEdge(rgw_current, heads.first, po));
      else {
        bool covered = false;
        for (const Node * nonroothead : nonrootheads)
          if (graph.hasEdge(rgw_current, nonroothead, po) &&
              graph.hasEdge(nonroothead, readnd, po)) {
            covered = true;
            break;
          }
        if (heads.first && graph.hasEdge(rgw_current, heads.first, po))
          covered = true;
        assert(covered);
      }
      #endif
      rgw_id++;
      continue;
    }
    assert(!graph.areOrdered(rgw_current, readnd, po));
    // Add the edge to readnd
    graph.addEdge(rgw_current, readnd, po);
    // This edge makes rgw_current the nonroot head
    #ifndef NDEBUG
    auto newheads = graph.getHeadWrites(readnd, po);
    assert(newheads.second.size() == 1 &&
           *(newheads.second.begin()) == rgw_current);
    assert(isGood(rgw_current, ann));
    #endif
    return {false, true}; // done, change
  }

  // Tried all and none worked
  assert(rgw_id == remoteGoodWrites.size());
  return {true, false}; // impossible
}

/* *************************** */
/* RULE 2                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> VCValClosure::ruleTwo
(const PartialOrder& po, const Node * readnd,
 const VCAnnotation::Ann& ann)
{
  assert(isRead(readnd));
  assert(readnd->getProcessID() == graph.starRoot() || ann.oneGW ||
         wBounds.count(readnd));
  // Rule2: there exists a write such that
  // (i) GOOD (ii) TAIL
  //
  // Note: if there is a local and a remote
  // GOOD write, make TAIL the one which
  // forces to add a weaker edge

  if (ann.oneGW && ann.loc != VCAnnotation::Loc::LOCAL) {
    /* ***************** */
    /* DATACENTRIC       */
    /* ***************** */
    // Rule2: r observes w ...
    // w -> w' implies r -> w'
    const Node *good = getGood(ann);
    assert(readnd->getProcessID() != graph.starRoot() &&
           "Rule2 done only for leaf read nodes");

    bool change = false;
    while (true) {
      bool changeThisIter = false;
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
        assert(!graph.areOrdered(readnd, roottail, po));
        graph.addEdge(readnd, roottail, po);
        change = true; changeThisIter = true;
      }
      for (const Node * nonroottail : nonroottails)
        if (nonroottail && graph.hasEdge(good, nonroottail, po)) {
          if (graph.hasEdge(nonroottail, readnd, po))
            return {true, false}; // impossible
          assert(!graph.hasEdge(readnd, nonroottail, po));
          graph.addEdge(readnd, nonroottail, po);
          change = true; changeThisIter = true;
        }
      if (!changeThisIter)
        return {true, false}; // impossible
    }

    return {false, change}; // done, change-bool
  }

  assert(!ann.oneGW || ann.loc == VCAnnotation::Loc::LOCAL);
  assert(readnd->getProcessID() != graph.starRoot() &&
         "Root doesn't need Rule2");

  /* ***************** */
  /* NONROOT   ANY     */
  /* ***************** */

  assert(ann.loc == VCAnnotation::Loc::ANY);
  assert(wBounds.count(readnd));
  #ifndef NDEBUG
  int& low = wBounds[readnd].first;
  #endif
  int& high = wBounds[readnd].second;
  const std::vector<const Node *>&
    wRemote = graph.wRoot.at(readnd->getEvent()->ml);
  assert(high < (int) wRemote.size());

  bool change = false;
  while (true) {
    auto tails = graph.getTailWrites(readnd, po);
    const Node *roottail = tails.first;
    if (roottail && isGood(roottail, ann)) {
      return {false, change}; // done, change-bool
    }
    // Since readnd is nonroot, all its nonroot tails
    // are ordered with (i.e. happen before) him
    const auto& nonroottails = tails.second;
    if (nonroottails.size() == 1 &&
        isGood(*(nonroottails.begin()), ann)) {
      return {false, change}; // done, change-bool
    }

    bool badNonrootTails = (nonroottails.size() > 1 ||
                            (nonroottails.size() == 1 &&
                             !isGood(*(nonroottails.begin()), ann)));
    #ifndef NDEBUG
    if (badNonrootTails)
      for (const Node * badNonrootTail : nonroottails)
        assert(!isGood(badNonrootTail, ann) &&
               graph.hasEdge(badNonrootTail, readnd, po));
    #endif

    // If nonroottails are bad they are lost
    // for good and we can focus only on root
    if (badNonrootTails) {
      if (wRemote.size() == 0 || low > high) {
        // Nonroottails are bad and there is no root tail to see
        return {true, false}; // impossible
      }
      assert(low <= high);
      assert(roottail && roottail == wRemote[high] &&
             !isGood(wRemote[high], ann));
      while (!isGood(wRemote[high], ann)) {
        if (low == high) {
          // Last root tail to see and it's bad
          return {true, false}; // impossible
        }
        assert(low < high && !graph.areOrdered(wRemote[high], readnd, po));
        high--;
      }
      assert(isGood(wRemote[high], ann) && (high + 1) < (int) wRemote.size() &&
             !isGood(wRemote[high + 1], ann) &&
             !graph.areOrdered(readnd, wRemote[high + 1], po));
      // We make wRemote[high] the good tail
      graph.addEdge(readnd, wRemote[high + 1], po);
      #ifndef NDEBUG
      auto newtails = graph.getTailWrites(readnd, po);
      assert(newtails.first == wRemote[high]);
      #endif
      return {false, true}; // done, change
    }
    // Nonroottails are not bad which means there are none currently,
    // there could be a good nonroot write possible to become tail,
    // (there definitely is some possible nonroot write - initial write),
    // so we have to check both options by another loop-iteration
    assert(nonroottails.size() == 0);
    assert(roottail && roottail == wRemote[high]);
    assert(!isGood(wRemote[high], ann) &&
           !graph.areOrdered(readnd, wRemote[high], po));
    graph.addEdge(readnd, wRemote[high], po);
    high--; // Original high now happens after readnd
    change = true;
  }
}

/* *************************** */
/* RULE 3                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> VCValClosure::ruleThree
(const PartialOrder& po, const Node * readnd,
 const VCAnnotation::Ann& ann)
{
  assert(isRead(readnd));
  assert(readnd->getProcessID() == graph.starRoot() || ann.oneGW ||
         wBounds.count(readnd));
  // Rule3: for every write such that
  // (i) BAD (ii) HEAD (iii) HB read
  // it HB some visible GOOD one
  //
  // Note: after Rule1, for every write
  // satisfying (i,ii,iii) and not HB
  // some visible GOOD write,
  // we need it to HB the 'other' tail
  // (which is GOOD after Rule2)

  if (ann.oneGW && ann.loc != VCAnnotation::Loc::LOCAL) {
    /* ***************** */
    /* DATACENTRIC       */
    /* ***************** */
    // Rule3: r observes w ...
    // w' -> r implies w' -> w
    const Node *good = getGood(ann);
    auto heads = graph.getHeadWrites(readnd, po);
    const auto& nonrootheads = heads.second;

    bool nonrootSituation = (nonrootheads.size() > 1 ||
                             (nonrootheads.size() == 1 &&
                              (*(nonrootheads.begin())) != good));

    if (!nonrootSituation)
      return {false, false}; // always done, no change

    #ifndef NDEBUG
    auto tails = graph.getTailWrites(readnd, po);
    assert(tails.first && tails.first == good);
    #endif

    bool change = false;
    for (const Node * badNonrootHead : nonrootheads) {
      // If any head doesn't have an edge
      // to the good root tail, add it
      if (graph.hasEdge(badNonrootHead, readnd, po) &&
          !graph.hasEdge(badNonrootHead, good, po)) {
        graph.addEdge(badNonrootHead, good, po);
        change = true;
      }
    }

    return {false, change}; // always done, change-bool

    // Root situation doesn't need to be dealt with
    // Max-min execution takes care of that no matter
    // whether the good write is root or nonroot
    // (cf Rule2 for leaf read, Rule1 for root read)
  }

  assert(!ann.oneGW || ann.loc == VCAnnotation::Loc::LOCAL);

  if (readnd->getProcessID() != graph.starRoot()) {
    /* ***************** */
    /* NONROOT   ANY     */
    /* ***************** */

    assert(ann.loc == VCAnnotation::Loc::ANY);
    auto heads = graph.getHeadWrites(readnd, po);
    const auto& nonrootheads = heads.second;

    bool nonrootSituation = (nonrootheads.size() > 1 ||
                             (nonrootheads.size() == 1 &&
                              !isGood(*(nonrootheads.begin()), ann)));

    if (!nonrootSituation)
      return {false, false}; // always done, no change

    assert(wBounds.count(readnd));
    int& high = wBounds[readnd].second;
    const std::vector<const Node *>&
      wRemote = graph.wRoot.at(readnd->getEvent()->ml);

    #ifndef NDEBUG
    assert(nonrootSituation);
    for (const Node * badNonrootHead : nonrootheads) {
      assert(!isGood(badNonrootHead, ann));
      assert(graph.hasEdge(badNonrootHead, readnd, po));
    }
    #endif

    // Get root tail, it is good (cf Rule2)
    const Node *roottail = wRemote[high];
    assert(isGood(roottail, ann));
    #ifndef NDEBUG
    auto tails = graph.getTailWrites(readnd, po);
    assert(tails.first && tails.first == roottail);
    #endif

    bool change = false;
    for (const Node * badNonrootHead : nonrootheads) {
      // If any head doesn't have an edge
      // to the good root tail, add it
      if (!graph.hasEdge(badNonrootHead, roottail, po)) {
        graph.addEdge(badNonrootHead, roottail, po);
        change = true;
      }
    }

    return {false, change}; // always done, change-bool

    // Root situation doesn't need to be dealt with
    // Max-min execution takes care of that no matter
    // whether a good tail (cf Rule2) is root or nonroot
  }

  assert(readnd->getProcessID() == graph.starRoot());

  if (ann.loc == VCAnnotation::Loc::LOCAL) {
    /* ***************** */
    /* ROOT   LOCAL      */
    /* ***************** */

    auto heads = graph.getHeadWrites(readnd, po);
    const auto& nonrootheads = heads.second;
    assert(heads.first && isGood(heads.first, ann));

    if (nonrootheads.size() == 0)
      return {false, false}; // always done, no change

    bool change = false;
    for (const Node * badNonrootHead : nonrootheads) {
      // Trivially bad since nonroot
      assert(!isGood(badNonrootHead, ann));
      // Not covered by the good root head
      assert(!graph.areOrdered(badNonrootHead, heads.first, po));
      // If HB readnd we have to make it
      // also HB the good root head
      if (graph.hasEdge(badNonrootHead, readnd, po)) {
        graph.addEdge(badNonrootHead, heads.first, po);
        change = true;
      }
    }

    // All bad heads are now covered, there may
    // be new nonroot heads now, but those can't
    // HB readnd otherwise they would cover
    // the previous nonroot heads
    #ifndef NDEBUG
    auto newheads = graph.getHeadWrites(readnd, po);
    for (const Node * badNonrootHead : newheads.second) {
      assert(!graph.areOrdered(badNonrootHead, readnd, po));
    }
    #endif

    return {false, change}; // always done, change-bool
  }

  assert(ann.loc == VCAnnotation::Loc::REMOTE);

  /* ***************** */
  /* ROOT   REMOTE     */
  /* ***************** */

  // This doesn't need to be dealt with
  // By Rule1 a visible good leaf write HB read
  // Max-min execution makes this write be the
  // observed one by the read

  #ifndef NDEBUG
  auto headsX = graph.getHeadWrites(readnd, po);
  const auto& nonrootheads = headsX.second;
  assert(nonrootheads.size() == 1 &&
         isGood(*(nonrootheads.begin()), ann) &&
         graph.hasEdge(*(nonrootheads.begin()), readnd, po));
  #endif

  return {false, false}; // always done, no change
}

/* *************************** */
/* RULES                       */
/* *************************** */

std::pair<bool, bool> VCValClosure::rules
(const PartialOrder& po, const Node * readnd,
 const VCAnnotation::Ann& ann)
{
  assert(readnd && isRead(readnd));
  //graph.to_dot(po,"");
  //readnd->dump();
  //ann.dump();
  // Update the bounds if it is a nonroot write
  if (readnd->getProcessID() != graph.starRoot() && !ann.oneGW)
    updateBounds(po, readnd);

  bool change = false;
  std::pair<bool, bool> res;
  // Rule1
  if (readnd->getProcessID() == graph.starRoot()) {
    // Leaf read doesn't need Rule1,
    // check paper for a proof
    res = ruleOne(po, readnd, ann);
    if (res.first) return {true, false};
    if (res.second) change = true;
  }
  // Rule2
  if (readnd->getProcessID() != graph.starRoot()) {
    // Root read doesn't need Rule2,
    // check paper for a proof
    res = ruleTwo(po, readnd, ann);
    if (res.first) return {true, false};
    if (res.second) change = true;
  }
  //Rule3
  res = ruleThree(po, readnd, ann);
  if (res.first) return {true, false};
  if (res.second) change = true;

  return {false, change};
}


/* *************************** */
/* VAL CLOSE                   */
/* *************************** */

void VCValClosure::valClose
(const PartialOrder& po, const Node * newread,
 const VCAnnotation::Ann * newann)
{
  prepare(po, newread, newann);

  bool change = false;
  if (newread) {
    auto res = rules(po, newread, *newann);
    if (res.first) { closed = false; return; }
    if (res.second) change = true;
  }
  if (newread && !change && graph.lessThanTwoLeavesWithRorW()) {
    // 1) this is not a closure after a lock-annotation
    // 2) there was no edge added to satisfy rules for newread
    // 3) we have one leaf with vis. event so no Maz orderings
    // Therefore we can state that the PO is already closed
    closed = true;
    return;
  }

  change = true;
  while (change) {
    change = false;
    for (const auto& key_ann : annotation) {
      assert(graph.closureSafeUntil.size() > key_ann.first.first);
      if (graph.closureSafeUntil[key_ann.first.first]
          < (int) key_ann.first.second) {
        auto res = rules(po, graph.getNode(key_ann.first), key_ann.second);
        if (res.first) { closed = false; return; }
        if (res.second) { change = true; }
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

void VCValClosure::valCloseLock(const PartialOrder& po,
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
