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

#include "ZAnnotation.h"


void ZAnnotation::add(const ZEvent *ev, const ZObs& obs)
{
  assert(isRead(ev));
  auto key = ZObs(ev->threadID(), ev->eventID());
  auto it = mapping.find(key);
  assert(it == mapping.end());
  mapping.emplace_hint(it, key, obs);
}


void ZAnnotation::add(const ZEvent *ev, const ZEvent *obsEv)
{
  assert(isRead(ev) && isWriteB(obsEv));
  add(ev, ZObs(obsEv->threadID(), obsEv->eventID()));
}


bool ZAnnotation::defines(unsigned thrid, unsigned evid) const
{
  auto key = ZObs(thrid, evid);
  return (mapping.find(key) != mapping.end());
}


bool ZAnnotation::defines(const ZEvent *ev) const
{
  assert(isRead(ev));
  return defines(ev->threadID(), ev->eventID());
}


const ZObs& ZAnnotation::getObs(unsigned thrid, unsigned evid) const
{
  auto key = ZObs(thrid, evid);
  auto it = mapping.find(key);
  assert(it != mapping.end());
  return it->second;
}


const ZObs& ZAnnotation::getObs(const ZEvent *ev) const
{
  assert(isRead(ev));
  return getObs(ev->threadID(), ev->eventID());
}


void ZAnnotation::setLastLock(const ZEvent *ev)
{
  assert(isLock(ev));
  auto it = lastlock.find(ev->ml);
  if (it != lastlock.end())
    it = lastlock.erase(it);
  lastlock.emplace_hint(it, ev->ml,
                        ZObs(ev->threadID(), ev->eventID()));
}


const ZObs& ZAnnotation::getLastLock(const ZEvent *ev) const
{
  assert(isLock(ev));
  auto it = lastlock.find(ev->ml);
  assert(it != lastlock.end());
  return it->second;
}


bool ZAnnotation::isLastLock(const ZEvent *ev) const
{
  assert(isLock(ev));
  auto it = lastlock.find(ev->ml);
  if (it == lastlock.end())
    return false;
  else
    return (it->second ==
            ZObs(ev->threadID(), ev->eventID()));
}


bool ZAnnotation::locationHasSomeLock(const ZEvent *ev) const
{
  assert(isLock(ev));
  return (lastlock.count(ev->ml));
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZAnnotation& annot) {
  out << "Annotation: {\n";
  for (auto& an : annot) {
    out << an.first.to_string() << "  observes::  ";
    out << an.second.to_string();
    out << "\n";
  }
  out << "}\n";
  return out;
}


void ZAnnotation::dump() const {
  llvm::errs() << *this;
}
