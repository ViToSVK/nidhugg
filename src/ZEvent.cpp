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
#include "ZHelpers.h"


ZEvent::ZEvent
(const IID<int> &iid, const CPid& cpid,
 int instruction_id, int event_id, int trace_id)
  : kind(Kind::DUMMY),
    _id(ZEventID(cpid, event_id)),
    _thread_id(-1), /*set at tracebuilder time*/
    _aux_id((cpid.is_auxiliary()) ? cpid.get_aux_index() : -1),
    _trace_id(trace_id),
    observed_trace_id(-1),
    write_other_trace_id(-1),
    write_other_ptr(nullptr), /*set at explorer.extend time*/
    childs_cpid(),
    fence(false),
    ml(SymAddr(SymMBlock::Stack(cpid.get_hash() % INT16_MAX, INT16_MAX), INT16_MAX), INT16_MAX),
    value(-1),
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
    _thread_id(-1),
    _aux_id(-1),
    _trace_id(-1),
    observed_trace_id(-1),
    write_other_trace_id(-1),
    write_other_ptr(nullptr),
    childs_cpid(),
    fence(false),
    ml(SymAddr(SymMBlock::Stack(INT16_MAX, INT16_MAX), INT16_MAX), INT16_MAX),
    value(-1),
    //
    iid(),
    size(-1),
    md(nullptr),
    may_conflict(false),
    instruction_id(-1)
{
  assert(initial);
}

// Returns a 'copy' of the event, with custom trace id
// The event will be part of replay_trace
ZEvent::ZEvent(const ZEvent& oth, int trace_id, bool keepvalue)
  : kind(oth.kind),
    _id(oth._id),
    _thread_id(oth._thread_id),
    _aux_id(oth._aux_id),
    _trace_id(trace_id),
    observed_trace_id(-1), /*not set for replay_trace events*/
    write_other_trace_id(-1), /*not set for replay_trace events*/
    write_other_ptr(nullptr), /*not set for replay_trace events*/
    childs_cpid(oth.childs_cpid),
    fence(oth.fence),
    ml(oth.ml),
    value(keepvalue ? oth.value : -1), /*not set for replay_trace events*/
    //
    iid(oth.iid), /*guide interpreter*/
    size(oth.size), /*guide interpreter*/
    md(nullptr),
    may_conflict(oth.may_conflict),
    instruction_id(oth.instruction_id) /*guide interpreter*/
{
  assert(!isInitial(oth) && "Should not copy initial event");
  assert(iid.get_pid() >= 0);
}


std::string ZEvent::to_string(bool write_cpid) const
{
  if (this->kind == ZEvent::Kind::INITIAL)
    return "-1::<>_[t-1,,a-1,,e-1] initial_event";

  std::stringstream res;

  res << trace_id() << "::";
  if (write_cpid)
    res << cpid() << "_";
  res << "[t" << thread_id() << "_a" << aux_id() << "_e" << event_id() << "]";
  switch(kind) {
   case ZEvent::Kind::DUMMY :
     res << " <s:" << size << ">";
     break;
   case ZEvent::Kind::READ :
     res << " read " << ml.addr.to_string()
         << " <O:" << observed_trace_id << "> <val:" << value << ">";
     break;
   case ZEvent::Kind::WRITEB :
     res << " writeB " << ml.addr.to_string() << " <val:" << value << ">";
     break;
   case ZEvent::Kind::WRITEM :
     res << " writeM " << ml.addr.to_string()
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
     res << " lock " << ml.addr.to_string()
         << " <O:" << observed_trace_id << ">";
     break;
   case ZEvent::Kind::M_UNLOCK :
     res << " unlock " << ml.addr.to_string();
     break;
   default :
     res << " unknown";
  }
  if (fence)
    res << " fence";

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


std::string trace_to_string(const std::vector<ZEvent>& trace)
{
  std::stringstream res;

  res << "TRACE::: " << trace.size() << " EVENTS\n";
  for (const auto& ev : trace)
    res << ev.to_string(true) << "\n";

  return res.str();
}


std::string trace_to_string(const std::vector<std::unique_ptr<ZEvent>>& trace)
{
  std::stringstream res;

  res << "TRACE::: " << trace.size() << " EVENTS\n";
  for (const auto& ev : trace) {
    assert(ev.get() && "No nullptrs stored in trace");
    res << ev->to_string(true) << "\n";
  }

  return res.str();
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const std::vector<ZEvent>& trace)
{
  out << trace_to_string(trace);
  return out;
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const std::vector<std::unique_ptr<ZEvent>>& trace)
{
  out << trace_to_string(trace);
  return out;
}


void dump_trace(const std::vector<ZEvent>& trace)
{
  llvm::errs() << trace << "\n";
}


void dump_trace(const std::vector<std::unique_ptr<ZEvent>>& trace)
{
  llvm::errs() << trace << "\n";
}
