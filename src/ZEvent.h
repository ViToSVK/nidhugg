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

#ifndef _Z_EVENT_H_
#define _Z_EVENT_H_

#if defined(HAVE_LLVM_IR_METADATA_H)
#include <llvm/IR/Metadata.h>
#elif defined(HAVE_LLVM_METADATA_H)
#include <llvm/Metadata.h>
#endif

#include <tuple>
#include <unordered_map>

#include "ZEventID.h"
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
  // Default constructor
  ZEvent(const IID<int> &iid, const CPid& cpid,
         int instruction_id, int event_id, int trace_id);
  // Constructor for initial event
  ZEvent(bool initial);

  // Returns a 'copy' of the event, with custom trace id
  // The event will be part of replay_trace
  ZEvent(const ZEvent& oth, int trace_id);

  ZEvent(const ZEvent& oth) = default;
  ZEvent(ZEvent&& oth) = default;
  ZEvent& operator=(const ZEvent& oth) = delete;
  ZEvent& operator=(ZEvent&& oth) = delete;

  enum class Kind {
    DUMMY, INITIAL,
    READ, WRITEB, WRITEM,
    SPAWN, JOIN,
    M_INIT, M_LOCK, M_UNLOCK, M_DESTROY
  } kind;

  /* A complex identifier of the thread that executed this event
   * + Sequential number (within the thread) of this event (starting with 0) */
 private:
  ZEventID _id;
 public:
  const ZEventID& id() const { return _id; }
  const CPid& cpid() const { return _id.cpid(); }
  int event_id() const { return _id.event_id(); }
  /* Thread ID in our partial order */
  mutable int _thread_id;
  int thread_id() const { return _thread_id; }
  /* Aux ID in our partial order */
  mutable int _aux_id;
  int aux_id() const { return _aux_id; }
  /* Trace ID of the event (index into the trace this event is part of) */
  int _trace_id;
  int trace_id() const { return _trace_id; }
  /* Trace ID of the observed event (lock observes unlock)
   * -1 means the initial event was observed */
  int observed_trace_id;
  /* Trace ID of its memory-write if this is a buffer-write
   * Trace ID of its buffer-write if this is a memory-write */
  mutable int write_other_trace_id;
  /* Point to this other event */
  mutable const ZEvent * write_other_ptr;
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
  int instruction_id;

  bool operator==(const ZEvent& c) const { return (id() == c.id()); }
  bool operator!=(const ZEvent &c) const { return (id() != c.id()); }
  bool operator<(const ZEvent &c) const { return (id() < c.id()); }
  bool operator<=(const ZEvent &c) const { return (id() <= c.id()); }
  bool operator>(const ZEvent &c) const { return (id() > c.id()); }
  bool operator>=(const ZEvent &c) const { return (id() >= c.id()); }

  std::string to_string(bool write_cpid = true) const;
  void dump() const;
};
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZEvent& ev);


class ZEventPtrComp {
 public:
  bool operator() (const ZEvent *e1, const ZEvent *e2) const {
    return (e1->operator<(*e2));
  }
};

std::string trace_to_string(const std::vector<ZEvent>& trace);
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const std::vector<ZEvent>& trace);
void dump_trace(const std::vector<ZEvent>& trace);

#endif // _Z_EVENT_H_
