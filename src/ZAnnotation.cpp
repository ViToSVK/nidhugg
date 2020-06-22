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


/* *************************** */
/* READS                       */
/* *************************** */

void ZAnnotation::add(const ZEventID& ev_id, const ZEventID& obs_id)
{
  assert(ev_id.event_id() >= 0 && !ev_id.cpid().is_auxiliary());
  auto it = read_mapping.find(ev_id);
  assert(it == read_mapping.end());
  read_mapping.emplace_hint(it, ev_id, obs_id);
}


void ZAnnotation::add(const ZEvent *ev, const ZEvent *obs)
{
  assert(isRead(ev) && isWriteB(obs));
  add(ev->id(), obs->id());
}


bool ZAnnotation::defines(const ZEventID& ev_id) const
{
  assert(ev_id.event_id() >= 0 && !ev_id.cpid().is_auxiliary());
  return (read_mapping.find(ev_id) != read_mapping.end());
}


bool ZAnnotation::defines(const ZEvent *ev) const
{
  bool res = defines(ev->id());
  assert(isRead(ev) || !res);
  return res;
}


const ZEventID& ZAnnotation::obs(const ZEventID& ev_id) const
{
  assert(ev_id.event_id() >= 0 && !ev_id.cpid().is_auxiliary());
  auto it = read_mapping.find(ev_id);
  assert(it != read_mapping.end());
  return it->second;
}


const ZEventID& ZAnnotation::obs(const ZEvent *ev) const
{
  assert(isRead(ev));
  return obs(ev->id());
}

/* *************************** */
/* LOCKS                       */
/* *************************** */

void ZAnnotation::lock_add(const ZEventID& ev_id, const ZEventID& obs_id)
{
  assert(ev_id.event_id() >= 0 && !ev_id.cpid().is_auxiliary());
  auto it = lock_mapping.find(ev_id);
  assert(it == lock_mapping.end());
  lock_mapping.emplace_hint(it, ev_id, obs_id);
}


void ZAnnotation::lock_add(const ZEvent *ev, const ZEvent *obs)
{
  assert(isLock(ev) && isUnlock(obs));
  lock_add(ev->id(), obs->id());
}


bool ZAnnotation::lock_defines(const ZEventID& ev_id) const
{
  assert(ev_id.event_id() >= 0 && !ev_id.cpid().is_auxiliary());
  return (lock_mapping.find(ev_id) != lock_mapping.end());
}


bool ZAnnotation::lock_defines(const ZEvent *ev) const
{
  bool res = lock_defines(ev->id());
  assert(isLock(ev) || !res);
  return res;
}


const ZEventID& ZAnnotation::lock_obs(const ZEventID& ev_id) const
{
  assert(ev_id.event_id() >= 0 && !ev_id.cpid().is_auxiliary());
  auto it = lock_mapping.find(ev_id);
  assert(it != lock_mapping.end());
  return it->second;
}


const ZEventID& ZAnnotation::lock_obs(const ZEvent *ev) const
{
  assert(isLock(ev));
  return lock_obs(ev->id());
}


int ZAnnotation::compare(const ZAnnotation &c) const
{
  if (read_size() < c.read_size()) { assert(to_string() != c.to_string()); return -1; }
  if (read_size() > c.read_size()) { assert(to_string() != c.to_string()); return 1; }
  if (lock_size() < c.lock_size()) { assert(to_string() != c.to_string()); return -1; }
  if (lock_size() > c.lock_size()) { assert(to_string() != c.to_string()); return 1; }
  assert(size() == c.size() && read_size() == c.read_size() && lock_size() == c.lock_size());

  // Reads
  auto cit = c.read_begin();
  for (auto it = read_begin(); it != read_end(); ++it) {
    assert(cit != c.read_end());
    // Same read
    const ZEventID& id = it->first;
    const ZEventID& cid = cit->first;
    int comp = id.compare(cid);
    if (comp != 0) { assert(to_string() != c.to_string()); return comp; }
    // Same observation-buffer-write
    const ZEventID& wid = it->second;
    const ZEventID& cwid = cit->second;
    comp = wid.compare(cwid);
    if (comp != 0) { assert(to_string() != c.to_string()); return comp; }
    //
    ++cit;
  }
  assert(cit == c.read_end());

  // Locks
  cit = c.lock_begin();
  for (auto it = lock_begin(); it != lock_end(); ++it) {
    assert(cit != c.lock_end());
    // Same lock
    const ZEventID& id = it->first;
    const ZEventID& cid = cit->first;
    int comp = id.compare(cid);
    if (comp != 0) { assert(to_string() != c.to_string()); return comp; }
    // Same observation-unlock
    const ZEventID& uid = it->second;
    const ZEventID& cuid = cit->second;
    comp = uid.compare(cuid);
    if (comp != 0) { assert(to_string() != c.to_string()); return comp; }
    //
    ++cit;
  }
  assert(cit == c.lock_end());

  assert(to_string() == c.to_string());
  return 0;
}


std::string ZAnnotation::to_string() const
{
  std::stringstream res;

  res << "Annotation-reads: {\n";
  for (const auto& an : *this) {
    res << an.first.to_string() << "  observes::  ";
    res << an.second.to_string();
    res << "\n";
  }
  res << "}\n";
  if (!lock_mapping.empty()) {
    res << "Annotation-locks: {\n";
    for (const auto& an : lock_mapping) {
      res << an.first.to_string() << "  observes::  ";
      res << an.second.to_string();
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
