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


bool ZAnnotationNeg::forbids_initial(const ZEvent *ev) const
{
  assert(is_read(ev) || is_lock(ev));
  return mapping.count(ev->id());
}


bool ZAnnotationNeg::forbids(const ZEvent *ev, const ZEvent *obs) const
{
  assert(obs && !is_initial(obs) &&
         "Call the special function for init event");
  assert((is_read(ev) && is_write(obs)) ||
         (is_lock(ev) && is_unlock(obs)));
  auto it = mapping.find(ev->id());
  if (it == mapping.end())
    return false;
  auto it2 = it->second.find(obs->cpid());
  if (it2 == it->second.end())
    return false;
  return (it2->second >= obs->event_id());
}


void ZAnnotationNeg::update(const ZEvent *ev, NegativeT&& upd)
{
  assert(is_read(ev) || is_lock(ev));
  auto it = mapping.find(ev->id());
  if (it == mapping.end()) {
    mapping.emplace_hint(it, ev->id(), upd);
    return;
  }
  assert(it != mapping.end());
  assert(it->second.size() <= upd.size());
  for (const auto& cpid_limit : upd) {
    assert(it->second.find(cpid_limit.first) == it->second.end() ||
           it->second.at(cpid_limit.first) <= cpid_limit.second);
    it->second[cpid_limit.first] = cpid_limit.second;
    assert(it->second.at(cpid_limit.first) == cpid_limit.second);
  }
}


std::string ZAnnotationNeg::to_string() const
{
  std::stringstream res;

  res << "Negative: {\n";
  for (auto& an : *this) {
    res << an.first.to_string() << "  forbids::  ";
    for (const auto& cpid_limit : an.second)
      res << cpid_limit.first.to_string() << "_"
          << cpid_limit.second << " ";
    res << "\n";
  }
  res << "}\n";
  return res.str();
}



void ZAnnotationNeg::dump() const {
  llvm::errs() << to_string();
}
