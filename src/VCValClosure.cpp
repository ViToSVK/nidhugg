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
                  (0, graph.processes.at(graph.starRoot()).size()-1));
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
  const std::vector<Node *>&
    wRemote = graph.processes.at(graph.starRoot());

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

    // Get heads, if |nonr| = 1 and good (triv. HB readnd), done
    // Otherwise nonr 'lost for good' (if |nonr| = 0, all covered)
    // Need root, in there using bounds, find first root write that is
    // good and visible (not covered by any nonroot head), and
    // add edge from this write to readnd, update bounds

  }

  assert(readnd->getProcessID() == graph.starRoot());

  if (ann.loc == VCAnnotation::Loc::LOCAL) {
    /* ***************** */
    /* ROOT   LOCAL      */
    /* ***************** */

    // Check root head, that's the only way to satisfy
  }

  assert(ann.loc == VCAnnotation::Loc::REMOTE);

  /* ***************** */
  /* ROOT   REMOTE     */
  /* ***************** */

  // Remote good writes are fully ordered, so quicksort them
  // and target either the head you found on
  // the first try (if good), or if nonr heads are bad, the
  // first one that is after any (which implies after all) of them

  // (*) Get heads, if |nonr| = 1 and good and edge to readnd, done
  // Otherwise, take the target in the sorted remote good writes
  // which has no edge to readnd, add this edge,
  // assert (get heads and check that this write is the single nonroot head)

  return {false, false};
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
    // Since readnd is nonroot, all its nonroot tails
    // are ordered with (i.e. happen before) him

    // Nonroot tails: if |.| > 1 assert all bad
    // if good done, if bad lost for good,
    // if |.| = 0 there is a chance (good one
    // above but not tail and could become)

    // (*) Get tails
    // if root not good, nonroot not good,
    // add edge readnd -> wbounds[high]
    // (assert it was the root tail),
    // decrease high, (*)

  }

  assert(readnd->getProcessID() == graph.starRoot());

  if (ann.loc == VCAnnotation::Loc::LOCAL) {
    /* ***************** */
    /* ROOT   LOCAL      */
    /* ***************** */

    // (*) Check root tail, after 1) the write above is good
    // if no root tail found, the write above is not tail,
    // so add edge from readnd to all nonroot tails, (*)
  }

  assert(ann.loc == VCAnnotation::Loc::REMOTE);

  /* ***************** */
  /* ROOT   REMOTE     */
  /* ***************** */

  // (*) Check nonroot tail(s), (if more assert all bad)
  // if bad, add edge from readnd to all of them, (*)

  return {false, false};
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

    // Get heads,
    // (because of Rule1 at most one of I, II can happen

    // (I) if bad root head and HB read: check root writes
    // (using wbounds) if there is some good, if there
    // isn't, get tails, assert nonroot good and HB readnd,
    // badroothead -> goodnonrtail, ++low (got covered)

    // (II) if bad nonroot head(s): assert all hb read,
    // wbounds[high] assert good, assert it's root tail
    // if any heads don't have, add -> goodroottail

  }

  assert(readnd->getProcessID() == graph.starRoot());

  if (ann.loc == VCAnnotation::Loc::LOCAL) {
    /* ***************** */
    /* ROOT   LOCAL      */
    /* ***************** */

    // (*) Get heads, check nonroot heads
    // if bad and HB readnd, because LOCAL they have
    // to be made to HB the root head
    // (which is also tail, assert good),
    // that covers them and there may be
    // new nonroot heads, but those can't
    // HB readnd otherwise they would cover
    // the previous nonroot heads, so done
  }

  assert(ann.loc == VCAnnotation::Loc::REMOTE);

  /* ***************** */
  /* ROOT   REMOTE     */
  /* ***************** */

  // Get heads, assert |nonr| = 1, good, hb read,
  // if root head exists, trivially bad,
  // get tails, assert |nonr| = 1, good,
  // if not already, rootbadhead -> nonrootgoodtail

  return {false, false};
}

/* *************************** */
/* RULES                       */
/* *************************** */

std::pair<bool, bool> VCValClosure::rules
(const PartialOrder& po, const Node * readnd,
 const VCAnnotation::Ann& ann)
{
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
