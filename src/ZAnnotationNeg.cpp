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

#include "ZAnnotationNeg.h"


bool ZAnnotationNeg::forbidsInitialEvent(const ZEvent *readev) const
{
  assert(isRead(readev) || isLock(readev));
  auto key = ZObs(readev->thread_id(), readev->event_id());
  return mapping.count(key);
}


bool ZAnnotationNeg::forbids(const ZEvent *readev, const ZEvent *writeev) const
{
  assert(writeev && !isInitial(writeev) &&
         "Call the special function for init event");
  assert((isRead(readev) && isWriteB(writeev)) ||
         (isLock(readev) && isUnlock(writeev)));
  assert(readev->aux_id() == -1 && writeev->aux_id() == -1);
  auto key = ZObs(readev->thread_id(), readev->event_id());
  auto it = mapping.find(key);
  if (it == mapping.end())
    return false;
  if (writeev->thread_id() >= it->second.size())
    return false;

  return (it->second[writeev->thread_id()]
          >= writeev->event_id());
}


void ZAnnotationNeg::update(const ZEvent *readev, std::vector<unsigned>&& newneg)
{
  assert(isRead(readev) || isLock(readev));
  auto key = ZObs(readev->thread_id(), readev->event_id());
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
  llvm::errs() << to_string();
}
