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

#include "ZEvent.h"


ZEvent::ZEvent
(const IID<int> &iid, const CPid& cpid,
 unsigned instruction_order, unsigned event_order, unsigned trace_id)
  : kind(Kind::DUMMY),
    cpid(cpid),
    childs_cpid(),
    fence(false),
    ml(SymAddr(SymMBlock::Stack(iid.get_pid(), 1337), 1337), 1337),
    value(-47),
    _thread_id(1337), /*set at PObuild time*/
    _aux_id((cpid.is_auxiliary()) ? cpid.get_aux_index() : -1),
    _event_id(event_order),
    _trace_id(trace_id),
    observed_trace_id(-1),
    write_other_trace_id(-1),
    write_other_ptr(nullptr), /*set at PObuild time*/
    //
    iid(iid),
    size(1),
    md(nullptr),
    may_conflict(false),
    instruction_order(instruction_order)
{
  assert(iid.get_pid() >= 0);
}


ZEvent::ZEvent(bool initial)
  : kind(Kind::INITIAL),
    ml(SymAddr(SymMBlock::Stack(INT16_MAX, INT16_MAX), INT16_MAX), INT16_MAX),
    value(-47),
    _thread_id(INT_MAX),
    _aux_id(-1),
    _event_id(INT_MAX),
    _trace_id(-1),
    observed_trace_id(-1),
    write_other_trace_id(-1),
    write_other_ptr(nullptr)
{
  assert(initial);
}


ZEvent::ZEvent(const ZEvent& oth, int trace_id, bool blank)
  : kind(oth.kind),
    cpid(oth.cpid),
    childs_cpid(oth.childs_cpid),
    fence(oth.fence),
    ml(oth.ml),
    value(-47),
    _thread_id(1337), /*set at PObuild time*/
    _aux_id(oth._aux_id),
    _event_id(oth._event_id),
    _trace_id(trace_id),
    observed_trace_id(-1),
    write_other_trace_id(-1),
    write_other_ptr(nullptr), /*set at PObuild time*/
    //
    iid(oth.iid), /*guide interpreter*/
    size(oth.size), /*guide interpreter*/
    md(nullptr),
    may_conflict(oth.may_conflict),
    instruction_order(oth.instruction_order) /*guide interpreter*/
{
  assert(iid.get_pid() >= 0);
  if (!blank) {
    // The observed trace ID stays the same
    // as we are reusing the whole trace
    observed_trace_id = oth.observed_trace_id;
    write_other_trace_id = oth.write_other_trace_id;
    value = oth.value;
  }
}


std::string ZEvent::to_string(bool write_cpid) const
{
  if (this->kind == ZEvent::Kind::INITIAL)
    return "-1_<>_[INT_MAX,-1,INT_MAX] initial_event";

  std::stringstream res;

  res << traceID() << "_";
  if (write_cpid)
    res << cpid << "_";
  res << "[" << threadID() << "," << auxID() << "," << eventID() << "]";
  switch(kind) {
   case ZEvent::Kind::DUMMY :
     res << " <s:" << size << ">";
     break;
   case ZEvent::Kind::READ :
     res << " read <- " << ml.addr.to_string()
         << " <O:" << observed_trace_id << "> <val:" << value << ">";
     break;
   case ZEvent::Kind::WRITEB :
     res << " writeB -> " << ml.addr.to_string() << " <val:" << value << ">";
     break;
   case ZEvent::Kind::WRITEM :
     res << " writeM -> " << ml.addr.to_string()
         << " <B:" << write_other_trace_id << "> <val:" << value << ">";
     break;
   case ZEvent::Kind::SPAWN :
     res << " spawn " << childs_cpid;
     break;
   case ZEvent::Kind::JOIN :
     res << " join " << childs_cpid;
     break;
   case ZEvent::Kind::M_INIT :
     res << " mutexinit " << ml.addr.to_string();
     break;
   case ZEvent::Kind::M_DESTROY :
     res << " mutexdestroy " << ml.addr.to_string();
     break;
   case ZEvent::Kind::M_LOCK :
     res << " lock " << ml.addr.to_string();
     break;
   case ZEvent::Kind::M_UNLOCK :
     res << " unlock " << ml.addr.to_string();
     break;
   default :
     res << " unknown";
  }

  return res.str();
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZEvent& ev)
{
  out << ev.to_string();
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
