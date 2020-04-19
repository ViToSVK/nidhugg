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

#include "ZEventID.h"


ZEventID::ZEventID(const CPid& cpid, int event_id)
  : _cpid(cpid),
    _event_id(event_id)
{
  _hash = compute_hash();
}


ZEventID::ZEventID(bool initial)
  : _cpid(),
    _event_id(-1)
{
  assert(initial);
  _hash = compute_hash();
}


std::size_t ZEventID::get_hash() const
{
  assert(compute_hash() == _hash);
  return _hash;
}


std::size_t ZEventID::compute_hash() const
{
  std::size_t res = get_cpid().get_hash();
  assert(res >= 10000 && res < 100000);
  res *= 1000;
  std::size_t ev = 999;
  assert(get_event_id() >= -1);
  if (_event_id >= 0)
    ev = (std::size_t) (_event_id % 1000);
  res += ev;
  assert(res >= 10000000 && res < 100000000);
  return res;
}


int ZEventID::compare(const ZEventID &c) const{
  std::size_t h = get_hash();
  std::size_t ch = c.get_hash();
  if(h < ch) return -1;
  if(h > ch) return 1;

  auto bad_hash = [&]()
  {
    llvm::errs() << "EventID hash clash: ";
    llvm::errs() << to_string() << " :: " << get_hash() << " --- ";
    llvm::errs() << c.to_string() << " :: " << c.get_hash() << "\n";
    assert(false && "EventID hash clash");
  };

  int cp = get_cpid().compare(c.get_cpid());
  if (cp < 0) { bad_hash(); return -1; }
  if (cp > 0) { bad_hash(); return 1; }

  if (get_event_id() < c.get_event_id()) { bad_hash(); return -1; }
  if (get_event_id() > c.get_event_id()) { bad_hash(); return 1; }

  return 0;
}


std::string ZEventID::to_string() const
{
  std::stringstream res;
  res << get_cpid().to_string() << "_";
  res << get_event_id();
  return res.str();
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZEventID& id)
{
  out << id.to_string();
  return out;
}


void ZEventID::dump() const
{
  llvm::errs() << *this << "\n";
}