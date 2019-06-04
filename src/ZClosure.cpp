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


std::pair<const ZEvent *, const ZEvent *>
ZClosure::getObs(const ZObs& obs)
{
  return {nullptr, nullptr};
}


/* *************************** */
/* RULE 1                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleOne
(const ZEvent *read, const ZObs& obs)
{
  return {false, false}; // done, no change
}


/* *************************** */
/* RULE 2                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleTwo
(const ZEvent *read, const ZObs& obs)
{
  return {false, false}; // done, no change
}


/* *************************** */
/* RULE 3                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleThree
(const ZEvent *read, const ZObs& obs)
{
  return {false, false}; // done, no change
}


/* *************************** */
/* RULES                       */
/* *************************** */

std::pair<bool, bool> ZClosure::rules
(const ZEvent *read, const ZObs& obs)
{
  assert(read && isRead(read));

  bool change = false;
  // Rule1 is done only the first time
  // Rule2
  auto res = ruleTwo(read, obs);
  if (res.first) return {true, false};
  if (res.second) change = true;
  //Rule3
  res = ruleThree(read, obs);
  if (res.first) return {true, false};
  if (res.second) change = true;

  return {false, change};
}


/* *************************** */
/* CLOSE                       */
/* *************************** */

bool ZClosure::close
(const ZEvent *newread)
{
  // Rules for new read
  if (newread) {
    auto res = rules(newread, an.getObs(newread));
    if (res.first) { return false; }
  }

  bool change = true;
  while (change) {
    change = false;
    for (const auto& read_obs : an) {
      auto read = gr.basis.getEvent(read_obs.first.thr, -1,
                                    read_obs.first.ev);
      if (!newread || newread != read) {
        auto res = rules(read, read_obs.second);
        if (res.first) { return false; }
        if (res.second) change = true;
      }
    }
    if (newread) {
      auto res = rules(newread, an.getObs(newread));
      if (res.first) { return false; }
      if (res.second) change = true;
    }
  }

  assert(!change);
  return true;
}


/* *************************** */
/* PRE-CLOSE                   */
/* *************************** */

void ZClosure::preClose
(const ZEvent *ev, const ZEvent *obsEv)
{
  assert(sameMl(ev, obsEv));
  assert((isRead(ev) && isWriteB(obsEv)) ||
         (isLock(ev) && isUnlock(obsEv)));

  if (isLock(ev)) {
    assert(!po.hasEdge(ev, obsEv));
    if (!po.hasEdge(obsEv, ev))
      po.addEdge(obsEv, ev);
    return;
  }

  assert(isRead(ev) && obsEv->write_other_ptr);
  if (ev->threadID() != obsEv->threadID()) {
    auto obsMem = obsEv->write_other_ptr;
    assert(isWriteM(obsMem));

    assert(!po.hasEdge(ev, obsMem));
    if (!po.hasEdge(obsMem, ev))
      po.addEdge(obsMem, ev);

    //TODO rest of rule1
  }
}


void ZClosure::preClose
(const ZEvent *ev, const ZObs& obs)
{
  preClose(ev, po.basis.getEvent(obs));
}
