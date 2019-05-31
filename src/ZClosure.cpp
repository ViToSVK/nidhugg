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

#include "ZClosure.h"


/* *************************** */
/* RULE 1                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleOne
(const ZEvent *read, const ZAnnotation::Ann& ann)
{
  return {false, false}; // done, no change
}


/* *************************** */
/* RULE 2                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleTwo
(const ZEvent *read, const ZAnnotation::Ann& ann)
{
  return {false, false}; // done, no change
}


/* *************************** */
/* RULE 3                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleThree
(const ZEvent *read, const ZAnnotation::Ann& ann)
{
  return {false, false}; // done, no change
}


/* *************************** */
/* RULES                       */
/* *************************** */

std::pair<bool, bool> ZClosure::rules
(const ZEvent *read, const ZAnnotation::Ann& ann)
{
  assert(read && isRead(read));

  bool change = false;
  // Rule1 is done only the first time
  // Rule2
  auto res = ruleTwo(read, ann);
  if (res.first) return {true, false};
  if (res.second) change = true;
  //Rule3
  res = ruleThree(read, ann);
  if (res.first) return {true, false};
  if (res.second) change = true;

  return {false, change};
}


/* *************************** */
/* CLOSE                       */
/* *************************** */

bool ZClosure::close
(const ZEvent *newread, const ZAnnotation::Ann *newobs)
{
  // Rule1 for new read + observation
  if (newread) {
    assert(newobs);
    auto res = ruleOne(newread, *newobs);
    if (res.first) { return false; }
  }

  bool change = true;
  while (change) {
    change = false;
    for (const auto& key_ann : an) {
      auto res = rules(gr.basis.getEvent(key_ann.first.first, -1, key_ann.first.second), key_ann.second);
      if (res.first) { return false; }
      if (res.second) change = true;
    }
    if (newread) {
      auto res = rules(newread, *newobs);
      if (res.first) { return false; }
      if (res.second) change = true;
    }
  }

  assert(!change);
  return true;
}


/* *************************** */
/* CLOSE LOCK                  */
/* *************************** */

bool ZClosure::closeLock
(const ZEvent *newlock, const ZEvent *lastunlock)
{
  assert(isLock(newlock) && isUnlock(lastunlock) &&
         sameMl(newlock, lastunlock));

  if (po.hasEdge(newlock, lastunlock)) {
    return false;
  } else if (po.hasEdge(lastunlock, newlock)) {
    return true;
  }

  po.addEdge(lastunlock, newlock);

  return close(nullptr, nullptr);
}

// test