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


int ZAnn::compare(const ZAnn& c) const
{
  if (value < c.value)
    return -1;
  else if (value > c.value)
    return 1;
  assert(false && "Comparing ZAnns with same values should not happen");
  return 0;
}


std::string ZAnn::to_string() const
{
  std::stringstream res;
  res << value;
  assert(!goodwrites.empty());
  for (const ZEventID& i : goodwrites)
    res << ":" << i.to_string();
  return res.str();
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZAnn& ann)
{
  out << ann.to_string();
  return out;
}


void ZAnn::dump() const {
  llvm::errs() << *this << "\n";
}


void ZAnnotation::add(const ZEventID& ev_id, const ZAnn& ann)
{
  assert(ev_id.event_id() >= 0 && !ev_id.cpid().is_auxiliary());
  auto it = mapping.find(ev_id);
  assert(it == mapping.end());
  mapping.emplace_hint(it, ev_id, ann);
}


bool ZAnnotation::defines(const ZEventID& ev_id) const
{
  assert(ev_id.event_id() >= 0 && !ev_id.cpid().is_auxiliary());
  return (mapping.find(ev_id) != mapping.end());
}


bool ZAnnotation::defines(const ZEvent *ev) const
{
  assert(is_read(ev));
  return defines(ev->id());
}


const ZAnn& ZAnnotation::ann(const ZEventID& ev_id) const
{
  assert(ev_id.event_id() >= 0 && !ev_id.cpid().is_auxiliary());
  auto it = mapping.find(ev_id);
  assert(it != mapping.end());
  return it->second;
}


const ZAnn& ZAnnotation::ann(const ZEvent *ev) const
{
  assert(is_read(ev));
  return ann(ev->id());
}


void ZAnnotation::set_last_lock(const ZEvent *ev)
{
  assert(is_lock(ev));
  auto it = lastlock.find(ev->ml());
  if (it != lastlock.end())
    it = lastlock.erase(it);
  lastlock.emplace_hint(it, ev->ml(), ev->id());
}


const ZEventID& ZAnnotation::last_lock(const ZEvent *ev) const
{
  assert(is_lock(ev));
  auto it = lastlock.find(ev->ml());
  assert(it != lastlock.end());
  return it->second;
}


bool ZAnnotation::is_last_lock(const ZEvent *ev) const
{
  assert(is_lock(ev));
  auto it = lastlock.find(ev->ml());
  if (it == lastlock.end())
    return false;
  else
    return (it->second == ev->id());
}


bool ZAnnotation::location_has_some_lock(const ZEvent *ev) const
{
  assert(is_lock(ev));
  return (lastlock.count(ev->ml()));
}


std::string ZAnnotation::to_string() const
{
  std::stringstream res;
  res << "Annotation: {\n";
  for (auto& an : *this) {
    res << an.first.to_string() << " ::obs:: ";
    res << an.second.to_string() << "\n";
  }
  res << "}\n";
  if (!lastlock.empty()) {
    res << "LastLock: {\n";
    for (const auto& last : lastlock) {
      res << last.first.addr.to_string() << " ::lockedby:: ";
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
  llvm::errs() << *this << "\n";
}
