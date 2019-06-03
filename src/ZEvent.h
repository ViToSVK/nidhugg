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

#ifndef _Z_EVENT_H_
#define _Z_EVENT_H_

#if defined(HAVE_LLVM_IR_METADATA_H)
#include <llvm/IR/Metadata.h>
#elif defined(HAVE_LLVM_METADATA_H)
#include <llvm/Metadata.h>
#endif

#include <tuple>
#include <unordered_map>

#include "CPid.h"
#include "TSOPSOTraceBuilder.h"
#include "SymEv.h"


/* Information about a (short) sequence of consecutive events by the
 * same thread. At most one event in the sequence may have conflicts
 * with other events, and if the sequence has a conflicting event,
 * it must be the LAST event in the sequence.
 */
class ZEvent {

 public:
  ZEvent() = delete;
  ZEvent(const IID<int> &iid, const CPid& cpid,
         unsigned instruction_order, unsigned event_order, unsigned trace_id)
   : kind(Kind::DUMMY),
    cpid(cpid),
    childs_cpid(),
    fence(false),
    ml(SymAddr(SymMBlock::Stack(iid.get_pid(), 1337), 1337), 1337),
    value(-47),
    _thread_id(1337), /*set at PObuild time*/
    _aux_id(cpid.get_aux_index()),
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
    instruction_order(instruction_order),
    aux_invisible()
    {
      assert(iid.get_pid() >= 0);
    }

  ZEvent(bool initial)
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

 private:
  // Returns a 'copy' of the event
  // If blank: the event will be a part of replay_trace
  // If not blank: the event is used as if it came from an interpreter
  ZEvent(const ZEvent& oth, int trace_id, bool blank)
    : kind(oth.kind),
    cpid(oth.cpid),
    childs_cpid(oth.childs_cpid),
    fence(oth.fence),
    ml(oth.ml),
    value(oth.value),
    _thread_id(1337), /*set at PObuild time*/
    _aux_id(oth.cpid.get_aux_index()),
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
    instruction_order(oth.instruction_order), /*guide interpreter*/
    aux_invisible(oth.aux_invisible) /*guide interpreter*/
    {
      assert(iid.get_pid() >= 0);
      if (!blank) {
        // The observed trace ID stays the same
        // as we are reusing the whole trace
        observed_trace_id = oth.observed_trace_id;
        write_other_trace_id = oth.write_other_trace_id;
      }
    }

 public:
  enum class Kind {
    DUMMY, INITIAL,
    READ, WRITEB, WRITEM,
    SPAWN, JOIN,
    M_INIT, M_LOCK, M_UNLOCK, M_DESTROY
  } kind;

  /* A complex identifier of the thread that executed this event */
  CPid cpid;
  /* CPid of the process that either:
   * 1) this event spawned (then this event is pthread_create) or
   * 2) this event joined (then this event is pthread_join) */
  CPid childs_cpid;
  /* Whether a fence happens immediately before the last
     (i.e. the visible) instruction of this event */
  bool fence;
  /* Memory location (if any) modified/read by this event */
  SymAddrSize ml;
  /* The value read/written by this event */
  int value;
  /* Thread ID in our partial order */
  mutable unsigned _thread_id;
  unsigned threadID() const { return _thread_id; }
  /* Aux ID in our partial order */
  int _aux_id;
  int auxID() const { return _aux_id; }
  /* Sequential number (within the thread) of this event,
   * also Event ID in our partial order
   * The first event of the thread is number 0 !!! */
  unsigned _event_id;
  unsigned eventID() const { return _event_id; }
  /* ID of the event (index into the trace this event is part of) */
  int _trace_id;
  int traceID() const { return _trace_id; }
  /* Trace ID of the observed event (lock observes unlock)
   * -1 means the initial event was observed */
  int observed_trace_id;
  /* Trace ID of its memory-write if this is a buffer-write
   * Trace ID of its buffer-write if this is a memory-write */
  mutable int write_other_trace_id;
  /* Point to this other event */
  mutable const ZEvent * write_other_ptr;

  ZEvent copy(int id, bool blank) const {
    return ZEvent(*this, id, blank);
  }

  /* Below for trace-builder */

  /* A simple identifier of the thread that executed this event */
  IID<int> iid;
  /* The number of instructions in this event. */
  int size;
  /* Metadata corresponding to the LAST instruction in this event. */
  const llvm::MDNode *md;
  /* Is it possible for the LAST instruction in this sequence to have a
   * conflict with an instruction in another event? */
  bool may_conflict;
  /* Sequential number (within the thread) of the LAST instruction in this event
   * The first instruction of the thread is number 1 !!! */
  unsigned instruction_order;
  /* Invisible instructions done on an auxiliary thread */
  std::unordered_map<int, unsigned> aux_invisible;

  bool operator==(const ZEvent& oth) const {
    return (threadID() == oth.threadID() &&
            auxID() == oth.auxID() &&
            eventID() == oth.eventID());
  }

  bool operator<(const ZEvent& oth) const {
    return std::tie(_thread_id, _aux_id, _event_id)
      < std::tie(oth._thread_id, oth._aux_id, oth._event_id);
  }

  std::string to_string(bool write_cpid) const;
  void dump() const;
};


class ZEventPtrComp {
 public:
  bool operator() (const ZEvent *e1, const ZEvent *e2) {
    return (e1->operator<(*e2));
  }
};


#endif // _Z_EVENT_H_
