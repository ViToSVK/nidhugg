/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2020 Viktor Toman
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

#include "ZAnnotationNeg.h"


bool ZAnnotationNeg::forbidsInitialEvent(const ZEvent *ev) const
{
  assert(isRead(ev) || isLock(ev));
  return mapping.count(ev->id());
}


bool ZAnnotationNeg::forbids(const ZEvent *ev, const ZEvent *obs) const
{
  assert(obs && !isInitial(obs) &&
         "Call the special function for init event");
  assert((isRead(ev) && isWriteB(obs)) ||
         (isLock(ev) && isUnlock(obs)));
  assert(ev->aux_id() == -1 && obs->aux_id() == -1);
  auto it = mapping.find(ev->id());
  if (it == mapping.end())
    return false;
  if (obs->thread_id() >= it->second.size())
    return false;

  return (it->second[obs->thread_id()]
          >= obs->event_id());
}


void ZAnnotationNeg::update(const ZEvent *ev, std::vector<unsigned>&& newneg)
{
  assert(isRead(ev) || isLock(ev));
  auto it = mapping.find(ev->id());
  if (it == mapping.end())
    mapping.emplace_hint(it, ev->id(), newneg);
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


std::string ZAnnotationNeg::to_string() const
{
  std::stringstream res;

  res << "Negative: {\n";
  for (auto& an : *this) {
    res << an.first.to_string() << "  forbids::  ";
    for (unsigned thr = 0; thr < an.second.size(); thr++)
      res << "[" << thr << ",-1," << an.second[thr] << "] ";
    res << "\n";
  }
  res << "}\n";
  return res.str();
}



void ZAnnotationNeg::dump() const {
  llvm::errs() << to_string() << "\n";
}
