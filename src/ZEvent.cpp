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

#include "ZEvent.h"


ZEvent::ZEvent
(const IID<int> &iid, const CPid& cpid,
 int instruction_id, int event_id, int trace_id)
  : kind(Kind::DUMMY),
    _id(ZEventID(cpid, event_id)),
    _trace_id(trace_id),
    _observed_trace_id(-1),
    _childs_cpid(),
    _ml(SymAddr(SymMBlock::Stack(cpid.get_hash() % INT16_MAX, INT16_MAX), INT16_MAX), INT16_MAX),
    _value(-1),
    //
    iid(iid),
    size(1),
    md(nullptr),
    may_conflict(false),
    instruction_id(instruction_id)
{
  assert(iid.get_pid() >= 0);
}


ZEvent::ZEvent(bool initial)
  : kind(Kind::INITIAL),
    _id(ZEventID(CPid(), -1)),
    _trace_id(-1),
    _observed_trace_id(-1),
    _childs_cpid(),
    _ml(SymAddr(SymMBlock::Stack(INT16_MAX, INT16_MAX), INT16_MAX), INT16_MAX),
    _value(0),
    //
    iid(),
    size(-1),
    md(nullptr),
    may_conflict(false),
    instruction_id(-1)
{
  assert(initial);
}


ZEvent::ZEvent(const ZEvent& oth, int trace_id,
               int observed_trace_id, int value)
  : kind(oth.kind),
    _id(oth.id()),
    _trace_id(trace_id),
    _observed_trace_id(observed_trace_id),
    _childs_cpid(oth.childs_cpid()),
    _ml(oth.ml()),
    _value(value),
    //
    iid(oth.iid), /*guide trace-builder*/
    size(oth.size), /*guide trace-builder*/
    md(nullptr),
    may_conflict(oth.may_conflict),
    instruction_id(oth.instruction_id) /*guide trace-builder*/
{
  assert(iid.get_pid() >= 0);
}


std::string ZEvent::to_string(bool write_cpid) const
{
  if (this->kind == ZEvent::Kind::INITIAL)
    return "-1::<>_-1 initial_event";

  std::stringstream res;

  res << trace_id() << "::";
  if (write_cpid)
    res << cpid() << "_";
  res << event_id();
  switch(kind) {
   case ZEvent::Kind::DUMMY :
     res << " <s:" << size << ">";
     break;
   case ZEvent::Kind::READ :
     res << " read " << ml().addr.to_string()
         << " <O:" << observed_trace_id() << "> <val:" << value() << ">";
     break;
   case ZEvent::Kind::WRITE :
     res << " write " << ml().addr.to_string() << " <val:" << value() << ">";
     break;
   case ZEvent::Kind::SPAWN :
     res << " spawn " << childs_cpid();
     break;
   case ZEvent::Kind::JOIN :
     res << " join " << childs_cpid();
     break;
   case ZEvent::Kind::M_INIT :
     res << " mutexinit " << ml().addr.to_string();
     break;
   case ZEvent::Kind::M_DESTROY :
     res << " mutexdestroy " << ml().addr.to_string();
     break;
   case ZEvent::Kind::M_LOCK :
     res << " lock " << ml().addr.to_string()
         << " <O:" << observed_trace_id() << ">";
     break;
   case ZEvent::Kind::M_UNLOCK :
     res << " unlock " << ml().addr.to_string();
     break;
   default :
     res << " unknown";
  }

  return res.str();
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZEvent& ev)
{
  out << ev.to_string(true);
  return out;
}


void ZEvent::dump() const
{
  llvm::errs() << *this << "\n";
}


void dumpTrace(const std::vector<ZEvent>& trace)
{
  llvm::errs() << "TRACE::: " << trace.size() << " EVENTS\n";
  for (const auto& ev : trace)
    ev.dump();
  llvm::errs() << "\n";
}
