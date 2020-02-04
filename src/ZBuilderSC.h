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

#ifndef __Z_BUILDER_SC_H__
#define __Z_BUILDER_SC_H__

#include <iostream>
#include <config.h>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <llvm/IR/Instructions.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/ExecutionEngine/GenericValue.h>

#include "Debug.h"
#include "DPORDriver.h"
#include "TSOTraceBuilder.h"
#include "Trace.h"
#include "ZEvent.h"


class ZBuilderSC : public TSOTraceBuilder {
  friend class ZExplorer;

  /* *************************** */
  /* SCHEDULING                  */
  /* *************************** */

  // We want to replay 'replay_trace'
  bool sch_replay;
  // We want to get an extension of our trace
  bool sch_extend;

  bool schedule_thread(int *proc, int *aux, unsigned p);
  bool schedule_arbitrarily(int *proc, int *aux);
  bool schedule_replay_trace(int *proc, int *aux);
  void update_prefix(unsigned p);

  void mayConflict(const SymAddrSize *ml = nullptr);

  /* *************************** */
  /* TRACES                      */
  /* *************************** */

  // The complete sequence of instructions executed since init of this TB
  // Format: vector of events; each event is a sequence of invisible instructions
  // followed by a single visible instruction (if any), all from the same thread
  std::vector<ZEvent> prefix;

  // We may have obtained this sequence as a constructor argument,
  // in such a case we first schedule in order to replay this entire sequence
  std::vector<ZEvent> replay_trace;

  // The index into prefix corresponding to the last event that was
  // scheduled. Has the value -1 when no events have been scheduled.
  // This is defined in TSOPSOTraceBuilder.h where we inherit from
  // I'm commenting it here so I'm aware of it
  // int prefix_idx = -1;
  // We might be receiving replay traces where a memory-write doesn't exactly
  // follow a buffer-write (eg b-write-ev7, other-thread-ev8, m-write-ev9)
  // We want to delay the buffer-write and play it right before the memory-write
  // (other-thread-ev7, b-write-ev8, m-write-ev9)
  int replay_trace_idx = -1;

  // Is there something in the extension to annotate?
  std::unordered_set<unsigned> somethingToAnnotate;
  // Does a thread end with a failed mutex lock attempt?
  // This can happen in case of a deadlock
  std::unordered_map<unsigned, SymAddrSize> endsWithLockFail;
  // Add lock event for every thread ending with a failed mutex lock attempt
  void add_failed_lock_attempts();
  // Trace index of the last write happening in the given location
  std::unordered_map<SymAddrSize, int> lastWrite;
  // Ipid -> Trace indices of writes in store queue of thread Ipid
  std::unordered_map<int, std::vector<int>> visibleStoreQueue;
  // During replay: in which ipids we are, for now, delaying a buffer-write
  std::unordered_set<unsigned> delayed_bwrite;


  ZEvent& curnode() {
    assert(0 <= prefix_idx);
    assert(prefix_idx < int(prefix.size()));
    return prefix[prefix_idx];
  };

  const ZEvent& curnode() const {
    assert(0 <= prefix_idx);
    assert(prefix_idx < int(prefix.size()));
    return prefix[prefix_idx];
  };

  /* *************************** */
  /* PARAMETER PASSING           */
  /* *************************** */
  unsigned star_root_index = 1;

 public:

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  // Use at the very beginning to get an initial trace
  ZBuilderSC(const Configuration &conf, llvm::Module *m,
             unsigned s_r_i);

  // Use when you want to do the following:
  // (step1) replay the trace tr
  // (step2) get a maximal extension
  ZBuilderSC(const Configuration &conf, llvm::Module *m,
             std::vector<ZEvent>&& tr);

  /* *************************** */
  /* CALLED FROM OUTSIDE         */
  /* *************************** */

  // Called from DPORDriver::run() at the beginning,
  // from here we create the explorer and explore
  virtual bool reset();

  // Called by Interpreter (Execution.cpp) in order to schedule one instruction
  // (not an entire event, just one (visible or invisible) instruction)
  virtual bool schedule(int *proc, int *aux, int *alt, bool *dryrun);

  // Called by Interpreter (Execution.cpp) to inform us
  // of a given special type of instruction happening
  virtual void refuse_schedule();  // called when scheduling got refused
  virtual void metadata(const llvm::MDNode *md);
  virtual IID<CPid> get_iid() const; // called from Interpreter::callFree
  //
  virtual void spawn();
  virtual void join(int tgt_proc);
  virtual void store(const SymData &sd, int val); // val
  virtual void store(const SymData &sd) {
    llvm::errs() << "Builder: Store should be reported along with the captured value\n";
    abort();
  }
  virtual void atomic_store(const SymData &sd);
  virtual void load(const SymAddrSize &ml, int val); // val
  virtual void load(const SymAddrSize &ml) {
    llvm::errs() << "Builder: Load should be reported along with the captured value\n";
    abort();
  }
  virtual void fence();
  virtual void mutex_init(const SymAddrSize &ml);
  virtual void mutex_destroy(const SymAddrSize &ml);
  virtual void mutex_lock(const SymAddrSize &ml);
  virtual void mutex_lock_fail(const SymAddrSize &ml);
  virtual void mutex_unlock(const SymAddrSize &ml);
  //
  virtual void mutex_trylock(const SymAddrSize &ml) {
    llvm::errs() << "Builder: No support for pthread_mutex_trylock\n";
    abort();
  }
  virtual bool cond_init(const SymAddrSize &ml) {
    llvm::errs() << "Builder: No support for cond_init\n";
    abort();
  }
  virtual bool cond_signal(const SymAddrSize &ml) {
    llvm::errs() << "Builder: No support for cond_signal\n";
    abort();
  }
  virtual bool cond_broadcast(const SymAddrSize &ml) {
    llvm::errs() << "Builder: No support for cond_broadcast\n";
    abort();
  }
  virtual bool cond_wait(const SymAddrSize &cond_ml, const SymAddrSize &mutex_ml) {
    llvm::errs() << "Builder: No support for cond_wait\n";
    abort();
  }
  virtual bool cond_awake(const SymAddrSize &cond_ml, const SymAddrSize &mutex_ml) {
    llvm::errs() << "Builder: No support for cond_awake\n";
    abort();
  }
  virtual int cond_destroy(const SymAddrSize &ml) {
    llvm::errs() << "Builder: No support for cond_destroy\n";
    abort();
  }
  virtual void compare_exchange
  (const SymData &sd, const SymData::block_type expected, bool success) {
    llvm::errs() << "Builder: No support for compare_exchange\n";
    abort();
  }
  virtual void full_memory_conflict() {
    llvm::errs() << "Builder: No support for full memory conflict: "
      << "(i) Do not assign to multiple variables in one command, "
      << "(ii) Do not use printf\n";
    abort();
  }

  // Called from ZExplorer on a TB created exclusively for this
  // Schedule entire replay_trace, then extend it, and return it
  std::pair<std::vector<ZEvent>, bool> extendGivenTrace();

  // We store an error trace (in their format) here
  // Trace *error_trace = nullptr;

  // Called by us when extending a trace, also called by DPORDriver
  // Inside the method, we translate an error trace (if we have any)
  // from our format into their format and then we return it
  virtual Trace *get_trace() const;
  virtual bool has_error() const { return (error_trace != nullptr)
                                            || (errors.size() > 0); };

};

#endif // __Z_BUILDER_SC_H__
