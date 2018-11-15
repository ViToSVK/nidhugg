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
(const PartialOrder& po, const Node * newread)
{
  if (newread) {
    // performing a new annotation
    // caches for all old nonroot annotations are set
    #ifndef NDEBUG
    for (auto& key_ann : annotation)
      if (!key_ann.second.ignore) {
        assert(graph.getNode(key_ann.first)->getProcessID() == graph.starRoot() ||
               wBounds.count(graph.getNode(key_ann.first)));
      }
    #endif
    if (newread->getProcessID() != graph.starRoot())
      prepareBounds(po, newread);
  } else {
    // either fully preparing an extension according to the
    // previous annotation, or closing after a lock annotation
    for (auto& key_ann : annotation)
      if (graph.getNode(key_ann.first)->getProcessID() != graph.starRoot())
        prepareBounds(po, graph.getNode(key_ann.first));
  }
}

void VCValClosure::prepareBounds
(const PartialOrder& po, const Node * readnd)
{
  assert(isRead(readnd->getEvent()));
  assert(readnd->getProcessID() != graph.starRoot());

  wBounds.emplace(readnd, std::pair<int, int>
                  (0, graph.wRoot.at(readnd->getEvent()->ml).size()-1));
  updateBounds(po, readnd);
}

void VCValClosure::updateBounds
(const PartialOrder& po, const Node * readnd)
{
  assert(isRead(readnd->getEvent()));
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
  assert(isRead(readnd->getEvent()));
  assert(readnd->getProcessID() == graph.starRoot() ||
         wBounds.count(readnd));
  // Rule1: there exists a write such that
  // (i) GOOD (ii) HEAD (iii) HB read
  //
  // Note: a local write is always
  // either *HEAD* or *covered*

  if (readnd->getProcessID() != graph.starRoot()) {
    /* ***************** */
    /* NONROOT   ANY     */
    /* ***************** */

    assert(ann.loc == VCAnnotation::Loc::ANY);
    // Since readnd is nonroot, all its nonroot heads
    // are ordered with (i.e. happen before) him
    auto heads = graph.getHeadWrites(readnd, po);
    const auto& nonrootheads = heads.second;
    if (nonrootheads.size() == 1 &&
        isGood(*(nonrootheads.begin()), ann)) {
      assert(graph.hasEdge(*(nonrootheads.begin()), readnd, po));
      return {false, false}; // done, no change
    }

    // Nonroot 'lost for good', need root head
    const Node *roothead = heads.first;
    if (roothead && isGood(roothead, ann)) {
      // We have good root head, now we
      // just need to check the HB edge
      assert(!graph.hasEdge(readnd, roothead, po));
      if (graph.hasEdge(roothead, readnd, po)) {
        // Edge is already there, we are done
        return {false, false}; // done, no change
      } else {
        // We have to add the edge
        graph.addEdge(roothead, readnd, po);
        // This edge doesn't change heads
        #ifndef NDEBUG
        auto newheads = graph.getHeadWrites(readnd, po);
        assert(newheads.first == roothead);
        assert(newheads.second.size() == heads.second.size());
        #endif
        return {false, true}; // done, change
      }
    }

    // Either no root head or a bad one, but some root write
    // still can be good, locate the first one like that
    assert(!roothead || !isGood(roothead, ann));
    assert(wBounds.count(readnd));
    int& low = wBounds[readnd].first;
    int& high = wBounds[readnd].second;
    const std::vector<const Node *>&
      wRemote = graph.wRoot.at(readnd->getEvent()->ml);
    if (wRemote.size() == 0)
      return {true, false}; // impossible
    assert(low <= high && low < (int) wRemote.size());
    if (!roothead && graph.hasEdge(wRemote[low], readnd, po)) {
      // wRemote[low] has edge but it's not head,
      // so it is covered, correct the bound
      low++;
      assert(low > high || !graph.areOrdered(wRemote[low], readnd, po));
    }
    if (roothead) {
      // Since there is a root head, it is bad
      // First make sure 'low' points to it
      if (roothead != wRemote[low]) {
        // Original 'low' is covered, correct the bound
        assert(graph.hasEdge(wRemote[low], readnd, po));
        low++;
        assert(low <= high && roothead == wRemote[low] &&
               !graph.areOrdered(roothead, readnd, po));
      }
      assert(roothead == wRemote[low] && !isGood(roothead, ann));
      // Increment 'low' to start search for a good write
      low++;
      assert(low > high || !graph.areOrdered(wRemote[low], readnd, po));
    }
    while (low <= high) {
      assert(!graph.areOrdered(wRemote[low], readnd, po));
      if (isGood(wRemote[low], ann)) {
        // Found first visible+good root write
        graph.addEdge(wRemote[low], readnd, po);
        // This edge makes wRemote[low] the root head
        // Nonroot heads could have got covered by wRemote[low]
        #ifndef NDEBUG
        auto newheads = graph.getHeadWrites(readnd, po);
        assert(newheads.first == wRemote[low]);
        #endif
        return {false, true}; // done, change
      } else
        low++;
    }

    // Didn't work with any root write
    return {true, false}; // impossible
  }

  assert(readnd->getProcessID() == graph.starRoot());

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
  // Remote good writes are fully ordered, quicksort them
  auto remoteGoodWrites = std::vector<const Node *>();
  auto itunord = graph.wNonrootUnord.find(readnd->getEvent()->ml);
  assert(itunord != graph.wNonrootUnord.end());
  for (auto& remoteWrite : itunord->second)
    if (isGood(remoteWrite, ann))
      remoteGoodWrites.push_back(remoteWrite);

  if (remoteGoodWrites.size() == 0)
    return {true, false}; // impossible
  #ifndef NDEBUG
  for (auto it1 = remoteGoodWrites.begin();
       it1 != remoteGoodWrites.end(); ++it1)
    for (auto it2 = remoteGoodWrites.begin();
         it2 != it1; ++it2)
      assert(graph.areOrdered(*it1, *it2, po));
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
      else
        for (const Node * nonroothead : nonrootheads)
          assert(graph.hasEdge(rgw_current, nonroothead, po));
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
  assert(isRead(readnd->getEvent()));
  assert(readnd->getProcessID() == graph.starRoot() ||
         wBounds.count(readnd));
  // Rule2: there exists a write such that
  // (i) GOOD (ii) TAIL
  //
  // Note: if there is a local and a remote
  // GOOD write, make TAIL the one which
  // forces to add a weaker edge

  if (readnd->getProcessID() != graph.starRoot()) {
    /* ***************** */
    /* NONROOT   ANY     */
    /* ***************** */

    assert(ann.loc == VCAnnotation::Loc::ANY);
    assert(wBounds.count(readnd));
    int& low = wBounds[readnd].first;
    int& high = wBounds[readnd].second;
    const std::vector<const Node *>&
      wRemote = graph.wRoot.at(readnd->getEvent()->ml);
    assert(high < (int) wRemote.size());

    bool change = false;
    while (true) {
      auto tails = graph.getTailWrites(readnd, po);
      const Node *roottail = tails.first;
      if (roottail && isGood(roottail, ann)) {
        return {false, change}; // always done, change-bool
      }
      // Since readnd is nonroot, all its nonroot tails
      // are ordered with (i.e. happen before) him
      const auto& nonroottails = tails.second;
      if (nonroottails.size() == 1 &&
          isGood(*(nonroottails.begin()), ann)) {
        return {false, change}; // always done, change-bool
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

      assert(low <= high && wRemote.size() > 0);
      // If nonroottails are bad they are lost
      // for good and we can focus only on root
      if (badNonrootTails) {
        assert(roottail && roottail == wRemote[high] &&
               !isGood(wRemote[high], ann));
        while (!isGood(wRemote[high], ann)) {
          assert(low < high);
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
        return {false, true}; // always done, change
      }
      // Nonroottails are not bad which means there are none currently,
      // there could be a good nonroot write possible to become tail,
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

  assert(readnd->getProcessID() == graph.starRoot());

  if (ann.loc == VCAnnotation::Loc::LOCAL) {
    /* ***************** */
    /* ROOT   LOCAL      */
    /* ***************** */

    bool change = false;
    while (true) {
      // After rule1, the write above is good
      auto tails = graph.getTailWrites(readnd, po);
      const Node *roottail = tails.first;
      if (roottail) {
        // The good write above is tail, done
        assert(isGood(roottail, ann));
        break;
      } else {
        // The write above is not tail,
        // add edge(s) from readnd to nonroot tail(s)
        assert(tails.second.size() > 0);
        for (const Node * nonroottail : tails.second) {
          assert(!graph.areOrdered(readnd, nonroottail, po));
          graph.addEdge(readnd, nonroottail, po);
        }
        change = true;
      }
    }
    return {false, change}; // always done, change-bool
  }

  assert(ann.loc == VCAnnotation::Loc::REMOTE);

  /* ***************** */
  /* ROOT   REMOTE     */
  /* ***************** */

  bool change = false;
  while (true) {
    // After rule 1, nonroot head is good
    // worst case it can be also tail
    auto tails = graph.getTailWrites(readnd, po);
    const auto& nonroottails = tails.second;
    assert(nonroottails.size() > 0);
    bool badtails = (nonroottails.size() > 1 ||
                     !isGood(*(nonroottails.begin()), ann));
    if (!badtails)
      break;
    // Bad tail(s), add edge from readnd to all of them
    for (const Node * badNonrootTail : nonroottails) {
      assert(!isGood(badNonrootTail, ann));
      assert(!graph.areOrdered(readnd, badNonrootTail, po));
      graph.addEdge(readnd, badNonrootTail, po);
    }
    change = true;
  }
  return {false, change}; // always done, change-bool
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
  assert(isRead(readnd->getEvent()));
  assert(readnd->getProcessID() == graph.starRoot() ||
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

  if (readnd->getProcessID() != graph.starRoot()) {
    /* ***************** */
    /* NONROOT   ANY     */
    /* ***************** */

    assert(ann.loc == VCAnnotation::Loc::ANY);
    auto heads = graph.getHeadWrites(readnd, po);
    const Node *roothead = heads.first;
    const auto& nonrootheads = heads.second;

    bool rootSituation = (roothead && !isGood(roothead, ann) &&
                          graph.hasEdge(roothead, readnd, po));
    bool nonrootSituation = (nonrootheads.size() > 1 ||
                             (nonrootheads.size() == 1 &&
                              !isGood(*(nonrootheads.begin()), ann)));

    assert(!rootSituation || !nonrootSituation); // After rule 1
    if (!rootSituation && !nonrootSituation)
      return {false, false}; // always done, no change

    assert(wBounds.count(readnd));
    int& low = wBounds[readnd].first;
    int& high = wBounds[readnd].second;
    const std::vector<const Node *>&
      wRemote = graph.wRoot.at(readnd->getEvent()->ml);

    if (rootSituation) {
      // First check root writes, if there is
      // any good, nothing needs to be done
      assert(roothead == wRemote[low] &&
             "Can't be low+1 because there's an edge");
      for (int idx = high; idx > low; --idx)
        if (isGood(wRemote[idx], ann))
          return {false, false}; // always done, no change

      // None was, to nonroot tail has to be good
      auto tails = graph.getTailWrites(readnd, po);
      const auto& nonroottails = tails.second;
      assert(nonroottails.size() == 1);
      const Node *singleNonrTail = *(nonroottails.begin());
      assert(isGood(singleNonrTail, ann) &&
             graph.hasEdge(singleNonrTail, readnd, po));
      assert(!graph.areOrdered(roothead, singleNonrTail, po));
      graph.addEdge(roothead, singleNonrTail, po);
      // Low got covered but we don't update the bounds in order
      // not to violate update assertions (low+1 has no edge)
      return {false, true}; // always done, change
    }

    assert(nonrootSituation);
    for (const Node * badNonrootHead : nonrootheads) {
      assert(!isGood(badNonrootHead, ann));
      assert(graph.hasEdge(badNonrootHead, readnd, po));
    }

    // Get root tail, it is good
    const Node *roottail = wRemote[high];
    assert(isGood(roottail, ann));
    #ifndef NDEBUG
    auto tails = graph.getTailWrites(readnd, po);
    assert(tails.first && tails.first == roottail);
    #endif

    bool change = false;
    for (const Node * badNonrootHead : nonrootheads) {
      // If any heads dont have an edge
      // to the good root tail, add it
      if (!graph.hasEdge(badNonrootHead, roottail, po)) {
        graph.addEdge(badNonrootHead, roottail, po);
        change = true;
      }
    }

    return {false, change}; // always done, change-bool
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

  auto heads = graph.getHeadWrites(readnd, po);
  const auto& nonrootheads = heads.second;
  assert(nonrootheads.size() == 1 &&
         isGood(*(nonrootheads.begin()), ann) &&
         graph.hasEdge(*(nonrootheads.begin()), readnd, po));

  if (!heads.first)
    return {false, false}; // always done, no change

  // If root head exists, trivially bad
  assert(!isGood(heads.first, ann));

  auto tails = graph.getTailWrites(readnd, po);
  const auto& nonroottails = tails.second;
  assert(nonroottails.size() == 1);
  const Node *singleNonrTail = *(nonroottails.begin());
  assert(isGood(singleNonrTail, ann) &&
         !graph.hasEdge(singleNonrTail, heads.first, po));

  if (!graph.hasEdge(heads.first, singleNonrTail, po)) {
    // This edge is needed for rule 3
    graph.addEdge(heads.first, singleNonrTail, po);
    return {false, true}; // always done, change
  }

  return {false, false}; // always done, no change
}

/* *************************** */
/* RULES                       */
/* *************************** */

std::pair<bool, bool> VCValClosure::rules
(const PartialOrder& po, const Node * readnd,
 const VCAnnotation::Ann& ann)
{
  assert(readnd && isRead(readnd->getEvent()));
  //graph.to_dot(po,"");
  //readnd->dump();
  //ann.dump();
  // Update the bounds if it is a nonroot write
  if (readnd->getProcessID() != graph.starRoot())
    updateBounds(po, readnd);

  bool change = false;
  // Rule1
  auto res = ruleOne(po, readnd, ann);
  if (res.first) return {true, false};
  if (res.second) change = true;
  // Rule2
  res = ruleTwo(po, readnd, ann);
  assert(!res.first);
  if (res.second) change = true;
  //Rule3
  res = ruleThree(po, readnd, ann);
  assert(!res.first);
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
  prepare(po, newread);

  bool change = true;
  while (change) {
    change = false;
    for (const auto& key_ann : annotation)
      if (!key_ann.second.ignore) {
        auto res = rules(po, graph.getNode(key_ann.first), key_ann.second);
        if (res.first) { closed = false; return; }
        if (res.second) change = true;
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
  assert(isLock(locknode->getEvent()) &&
         isUnlock(lastunlocknode->getEvent()) &&
         locknode->getEvent()->ml == lastunlocknode->getEvent()->ml);

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
