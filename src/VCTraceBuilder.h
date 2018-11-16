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

#ifndef __VC_TRACE_BUILDER_H__
#define __VC_TRACE_BUILDER_H__

#include <iostream>
#include <config.h>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <llvm/IR/Instructions.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/ExecutionEngine/GenericValue.h>

#include "Debug.h"
#include "DPORDriver.h"
#include "TSOTraceBuilder.h"
#include "Trace.h"
#include "VCEvent.h"

class VCTraceBuilder : public TSOTraceBuilder {

 public:

  /* *************************** */
  /* INTERNALS                   */
  /* *************************** */

  const Configuration &config;
  llvm::Module *M;

 private:

  /* *************************** */
  /* SCHEDULING                  */
  /* *************************** */

  // We want to get some initial trace
  bool sch_initial;
  // We want to replay 'replay_trace'
  bool sch_replay;
  // We want to get an extension of our trace
  bool sch_extend;

  bool schedule_thread(int *proc, unsigned p);
  bool schedule_arbitrarily(int *proc);
  bool schedule_replay_trace(int *proc);
  void update_prefix(unsigned p);

  // Size of the map is always 0 or 1, when 1:
  // key: id of process in CS; value: how deep in CS
  std::unordered_map<int, int> in_critical_section;

  // This is the currently executed instruction
  const llvm::Instruction *current_inst = nullptr;
  void mayConflict(const SymAddrSize *ml = nullptr);

  /* *************************** */
  /* TRACES                      */
  /* *************************** */

  // We consider the global's initialization event
  // It 'writes' initial value 0 to all global variables
  // const VCEvent initial_event = VCEvent(IID<IPid>(0,0));

  // The complete sequence of instructions executed since init of this TB
  // Format: vector of events; each event is a sequence of invisible instructions
  // followed by a single visible instruction (if any), all from the same thread
  // (we don't save anything about the invisible ones, just how many there are)
  std::vector<VCEvent> prefix;

  // We may have obtained this sequence as a constructor argument,
  // in such a case we first schedule in order to replay this entire sequence
  std::vector<VCEvent> replay_trace;

  // The set of threads that have executed a new so far unannotated
  // (visible) read; we refuse to schedule any such thread any further, we want to
  // annotate the new read and only next time allow the thread to progress further
  std::unordered_set<int> threads_with_unannotated_read;

  // The index into prefix corresponding to the last event that was
  // scheduled. Has the value -1 when no events have been scheduled.
  // This is defined in TSOPSOTraceBuilder.h where we inherit from
  // I'm commenting it here so I'm aware of it
  // int prefix_idx = -1;

  // Number of executed instructions since init of this TB
  unsigned executed_instr = 0;

  VCEvent& curnode() {
    assert(0 <= prefix_idx);
    assert(prefix_idx < int(prefix.size()));
    return prefix[prefix_idx];
  };

  const VCEvent& curnode() const {
    assert(0 <= prefix_idx);
    assert(prefix_idx < int(prefix.size()));
    return prefix[prefix_idx];
  };

  /* *************************** */
  /* PARAMETER PASSING           */
  /* *************************** */
  unsigned star_root_index = 1;
  bool previous_mutation_process_first = true;
  bool root_before_nonroots = true;

 public:

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  // Use at the very beginning to get an initial trace
  VCTraceBuilder(const Configuration &conf, llvm::Module *m,
                 unsigned s_r_i, bool p_m_p_f, bool r_b_n)
  : TSOTraceBuilder(conf), config(conf), M(m),
    sch_initial(true), sch_replay(false), sch_extend(false),
    in_critical_section(),
    star_root_index(s_r_i),
    previous_mutation_process_first(p_m_p_f),
    root_before_nonroots(r_b_n)
    {}

  // Use when you want to do the following:
  // (step1) replay the trace tr
  // (step2) extend it until each thread has a new
  //         unannotated read (or we can not extend further)
  VCTraceBuilder(const Configuration &conf,
                 llvm::Module *m, std::vector<VCEvent>&& tr,
                 const std::unordered_set<int>& unannot)
  : TSOTraceBuilder(conf), config(conf), M(m),
    sch_initial(false), sch_replay(true), sch_extend(false),
    in_critical_section(),
    replay_trace(std::move(tr)),
    threads_with_unannotated_read(unannot)
    {
      prefix.reserve(replay_trace.size() + 16);
    }

  /* *************************** */
  /* CALLED FROM OUTSIDE         */
  /* *************************** */

  // Called from DPORDriver::run() at the beginning,
  // from here we create the explorer and explore
  virtual bool reset();

  // Called by Interpreter (Execution.cpp) while executing the scheduled instruction
  virtual void executing_instruction(const llvm::Instruction *Instr) {
    current_inst = Instr;
    executed_instr++;
    // Since the clock part is refactored, I'm commenting the following out for now
    // unsigned p = curnode().iid.get_pid();
    // curnode().order = threads[p].clock[p];
  }

  // Called by Interpreter (Execution.cpp) in order to schedule one instruction
  // (not an entire event, just one (visible or invisible) instruction)
  virtual bool schedule(int *proc, int *aux, int *alt, bool *dryrun);

  // Called by Interpreter (Execution.cpp) to inform us
  // of a given special type of instruction happening
  virtual void refuse_schedule();  // called when scheduling got refused
  virtual void metadata(const llvm::MDNode *md);
  virtual IID<CPid> get_iid() const; // called from Interpreter::callFree
  virtual void spawn();
  virtual void join(int tgt_proc);
  virtual void atomic_store(const SymData &sd, int val); // WRITE (val)
  virtual void load(const SymAddrSize &ml, int val); // READ (val)
  virtual void fence();
  virtual void mutex_init(const SymAddrSize &ml);
  virtual void mutex_destroy(const SymAddrSize &ml);
  virtual void mutex_unlock(const SymAddrSize &ml); // WRITE (val=ipid << 16 + event_order)
  virtual void mutex_lock(const SymAddrSize &ml); // READ (val)
  virtual void mutex_lock_fail(const SymAddrSize &ml);
  virtual void mutex_trylock(const SymAddrSize &ml) {
    llvm::errs() << "No support for pthread_mutex_trylock\n";
    abort();
  }
  virtual void full_memory_conflict() {
    llvm::errs() << "No support for full memory conflict\n";
    abort();
  }

  // Called from VCExplorer on a TB created exclusively for this
  // Schedule entire replay_trace, then extend it, and return it
  std::pair<std::vector<VCEvent>,
    std::unordered_map<int, int>> extendGivenTrace();

  // We store an error trace (in their format) here
  // Trace *error_trace = nullptr;

  // Called by us when extending a trace, also called by DPORDriver
  // Inside the method, we translate an error trace (if we have any)
  // from our format into their format and then we return it
  virtual Trace *get_trace() const;
  virtual bool has_error() const { return (error_trace != nullptr)
                                            || (errors.size() > 0); };

};

#endif
