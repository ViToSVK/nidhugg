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


class ZEventAtomicInfo {
 public:
  ZEventAtomicInfo();
  ZEventAtomicInfo(const ZEventAtomicInfo&) = default;
  ZEventAtomicInfo(ZEventAtomicInfo&&) = default;
  ZEventAtomicInfo& operator=(const ZEventAtomicInfo&) = delete;
  ZEventAtomicInfo& operator=(ZEventAtomicInfo&&) = delete;

  bool is_part_of_cas() const { return _is_part_of_cas; }
  bool is_part_of_rmw() const { return _is_part_of_rmw; }
  int cas_compare_val() const { return _cas_compare_val; }
  int cas_exchange_val() const { return _cas_exchange_val; }

  void set_read_of_cas(int compare_val, int exchange_val);
  void set_write_of_cas();
  void set_read_of_rmw();
  void set_write_of_rmw();

 private:
  bool _is_part_of_cas;
  bool _is_part_of_rmw;
  int _cas_compare_val;
  int _cas_exchange_val;
};

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

  // Copy of the event, with custom-set trace id, observed trace id, value
  ZEvent(const ZEvent& oth, int trace_id,
         int observed_trace_id, int value);

  ZEvent(const ZEvent& oth) = default;
  ZEvent(ZEvent&& oth) = default;
  ZEvent& operator=(const ZEvent& oth) = delete;
  ZEvent& operator=(ZEvent&& oth) = delete;

  enum class Kind {
    DUMMY, INITIAL,
    READ, WRITE,
    SPAWN, JOIN,
    M_INIT, M_LOCK, M_UNLOCK, M_DESTROY
  } kind;

  /* A complex identifier of the thread that executed this event
   * + Sequential number (within the thread) of this event (starting with 0) */
  ZEventID _id;
  const ZEventID& id() const { return _id; }
  const CPid& cpid() const { return _id.cpid(); }
  int event_id() const { return _id.event_id(); }
  /* Trace ID of the event (index into the trace this event is part of) */
  int _trace_id;
  int trace_id() const { return _trace_id; }
  /* Trace ID of the observed event (lock observes unlock)
   * -1 means the initial event was observed */
  int _observed_trace_id;
  int observed_trace_id() const { return _observed_trace_id; }
  /* CPid of the process that either:
   * 1) this event spawned (then this event is pthread_create) or
   * 2) this event joined (then this event is pthread_join) */
  CPid _childs_cpid;
  const CPid& childs_cpid() const { return _childs_cpid; }
  /* Memory location (if any) modified/read by this event */
  SymAddrSize _ml;
  const SymAddrSize& ml() const { return _ml; }
  /* The value read/written by this event */
  int _value;
  int value() const { return _value; }
  /* Info whether this event is part of CAS or RMW */
  ZEventAtomicInfo _atomic_info;
  bool is_read_of_cas() const;
  bool is_write_of_cas() const;
  bool is_read_of_rmw() const;
  bool is_write_of_rmw() const;
  bool is_read_of_atomic_event() const;
  bool is_write_of_atomic_event() const;
  void set_read_of_cas(int compare_val, int exchange_val);
  void set_write_of_cas();
  void set_read_of_rmw();
  void set_write_of_rmw();
  int cas_compare_val() const;
  int cas_exchange_val() const;

  /* Guide trace-builder */

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

  std::string to_string(bool write_cpid) const;
  void dump() const;
};


class ZEventPtrComp {
 public:
  bool operator() (const ZEvent *e1, const ZEvent *e2) const {
    return (e1->operator<(*e2));
  }
};
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZEvent& ev);


void dump_trace(const std::vector<ZEvent>& trace);

#endif // _Z_EVENT_H_
