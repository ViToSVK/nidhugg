/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2018 Viktor Toman
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

#ifndef _VC_EVENT_H_
#define _VC_EVENT_H_

#if defined(HAVE_LLVM_IR_METADATA_H)
#include <llvm/IR/Metadata.h>
#elif defined(HAVE_LLVM_METADATA_H)
#include <llvm/Metadata.h>
#endif

#include <llvm/IR/Instruction.h>
#include <llvm/IR/Instructions.h>

#include "CPid.h"
#include "TSOPSOTraceBuilder.h"
#include "SymEv.h"

/* An identifier for a thread. An index into this->threads.
 *
 * Even indexes are for real threads. Odd indexes i are for
 * auxiliary threads corresponding to the real thread at index i-1.
 */
//typedef int IPid;

/* Information about a (short) sequence of consecutive events by the
 * same thread. At most one event in the sequence may have conflicts
 * with other events, and if the sequence has a conflicting event,
 * it must be the first event in the sequence.
 */
class VCEvent {
 public:

  VCEvent(const IID<int> &iid, const CPid& cpid,
					unsigned instruction_order, unsigned event_order,
					unsigned id)
    : iid(iid), cpid(cpid), childs_cpid(),
	    size(1), md(0), instruction(0),
		  may_conflict(false),
		  ml(SymAddr(SymMBlock::Stack(iid.get_pid(), 47), 47), 4),
		  value(0),
		  instruction_order(instruction_order),
		  event_order(event_order),
		  id(id)
      { assert(iid.get_pid() >= 0); }

  VCEvent(const IID<int> &iid, const CPid& cpid)
		: VCEvent(iid, cpid, 0, -1, -1) {}
	
  VCEvent(const IID<int> &iid)
		: VCEvent(iid, CPid(), 0, -1, -1) {}

  /* A simple identifier of the thread that executed this event */
  IID<int> iid;
  /* A complex identifier of the thread that executed this event */
  CPid cpid;
  /* CPid of the process that either:
   * 1) this event spawned (then this event is pthread_create) or
   * 2) this event joined (then this event is pthread_join) */
  CPid childs_cpid;
  /* The number of instructions in this event. */
  int size;
  /* Metadata corresponding to the LAST instruction in this event. */
  const llvm::MDNode *md;
  /* LAST instruction (the only visible one if any) in this event */
  const llvm::Instruction *instruction;
  /* Is it possible for the last instruction in this sequence to have a
   * conflict with an instruction in another event? */
  bool may_conflict;
  /* Memory location (if any) modified/read by this event */
  SymAddrSize ml;
	/* Value (if any) stored/loaded by this event */
	int value;
	/* Sequential number (within the thread) of the LAST instruction in this event
	 * The first instruction of the thread is number 1 !!! */
	unsigned instruction_order;
	/* Sequential number (within the thread) of this event
	 * The first event of the thread is number 0 !!! */
	unsigned event_order;
  /* ID of the event (index into the trace this event is part of) */
  unsigned id;

  // return a copy of this event without any computed information
  // contained in it -- that is copy only the basic attributes (iid, size,...)
  // but do not copy the branches, alts, etc..
  VCEvent blank_copy() const {
    VCEvent tmp(iid);

		tmp.cpid = cpid;
		tmp.size = size;
		tmp.md = md;
		tmp.instruction = instruction;
		// orders stay the same since the basis
		// up until this event stays the same
		tmp.instruction_order = instruction_order;
		tmp.event_order = event_order;
    
    return tmp;
  }

  void dump() const;
	
};

#endif // _VC_EVENT_H_
