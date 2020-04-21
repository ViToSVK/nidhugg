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

#include "ZAnnotation.h"


void ZAnnotation::add(const ZEventID& ev_id, const ZEventID& obs_id)
{
  assert(ev_id.event_id() >= 0 && ev_id.cpid().get_aux_index() == -1);
  auto it = mapping.find(ev_id);
  assert(it == mapping.end());
  mapping.emplace_hint(it, ev_id, obs_id);
}


void ZAnnotation::add(const ZEvent *ev, const ZEvent *obs)
{
  assert(isRead(ev) && isWrite(obs));
  add(ev->id(), obs->id());
}


bool ZAnnotation::defines(const ZEventID& ev_id) const
{
  assert(ev_id.event_id() >= 0 && ev_id.cpid().get_aux_index() == -1);
  return (mapping.find(ev_id) != mapping.end());
}


bool ZAnnotation::defines(const ZEvent *ev) const
{
  assert(isRead(ev));
  return defines(ev->id());
}


const ZEventID& ZAnnotation::obs(const ZEventID& ev_id) const
{
  assert(ev_id.event_id() >= 0 && ev_id.cpid().get_aux_index() == -1);
  auto it = mapping.find(ev_id);
  assert(it != mapping.end());
  return it->second;
}


const ZEventID& ZAnnotation::obs(const ZEvent *ev) const
{
  assert(isRead(ev));
  return obs(ev->id());
}


void ZAnnotation::set_last_lock(const ZEvent *ev)
{
  assert(isLock(ev));
  auto it = lastlock.find(ev->ml());
  if (it != lastlock.end())
    it = lastlock.erase(it);
  lastlock.emplace_hint(it, ev->ml(), ev->id());
}


const ZEventID& ZAnnotation::last_lock(const ZEvent *ev) const
{
  assert(isLock(ev));
  auto it = lastlock.find(ev->ml());
  assert(it != lastlock.end());
  return it->second;
}


bool ZAnnotation::is_last_lock(const ZEvent *ev) const
{
  assert(isLock(ev));
  auto it = lastlock.find(ev->ml());
  if (it == lastlock.end())
    return false;
  else
    return (it->second == ev->id());
}


bool ZAnnotation::location_has_some_lock(const ZEvent *ev) const
{
  assert(isLock(ev));
  return (lastlock.count(ev->ml()));
}


std::string ZAnnotation::to_string() const
{
  std::stringstream res;

  res << "Annotation: {\n";
  for (auto& an : *this) {
    res << an.first.to_string() << "  observes::  ";
    res << an.second.to_string();
    res << "\n";
  }
  res << "}\n";
  if (!lastlock.empty()) {
    res << "LastLock: {\n";
    for (const auto& last : lastlock) {
      res << last.first.to_string() << "  is locked by::  ";
      res << last.second.to_string();
      res << "\n";
    }
    res << "}\n";
  }
  return res.str();
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZAnnotation& annot)
{
  out << annot.to_string();
  return out;
}


void ZAnnotation::dump() const {
  llvm::errs() << *this;
}
