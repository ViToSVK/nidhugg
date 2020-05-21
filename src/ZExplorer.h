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

#ifndef __Z_EXPLORER_H__
#define __Z_EXPLORER_H__

#include <list>
#include <memory>
#include <time.h>

#include "ZBuilderSC.h"
#include "ZTrace.h"
#include "ZClosure.h"
#include "ZLinearization.h"


class ZExplorer {

  TSOPSOTraceBuilder * original_TB = nullptr;

  ZTrace * initial;

  bool info = false;

  class TraceExtension {
   public:
    TraceExtension() = default;
    TraceExtension(std::vector<ZEvent>&& extension,
                   bool some_to_ann, bool assume_blocked);

    TraceExtension(std::pair<std::vector<ZEvent>&&, bool>&& ext_sometoann);

    bool empty() const { return (trace.empty()); }

    std::vector<ZEvent> trace;
    bool something_to_annotate;
    bool has_assume_blocked_thread;
    bool has_error;
  };


  /* *************************** */
  /* ALGORITHM                   */
  /* *************************** */

 public:

  bool explore();

  void print_stats() const;

 private:

  bool explore_rec(ZTrace& ann_trace);

  bool mutate_read(const ZTrace& ann_trace, const ZEvent *read);

  bool mutate_lock(const ZTrace& ann_trace, const ZEvent *lock);

  bool close_po
  (const ZTrace& ann_trace, const ZEvent *read_lock,
   ZAnnotation&& mutated_annotation, ZPartialOrder&& mutated_po,
   bool mutation_follows_current_trace);

  bool extend_and_recur
  (const ZTrace& parent_trace, ZAnnotation&& mutated_annotation,
   ZPartialOrder&& mutated_po, bool mutation_follows_current_trace);

  TraceExtension reuse_trace
  (const ZTrace& parent_trace, const ZAnnotation& mutated_annotation);

  TraceExtension extend_trace(std::vector<ZEvent>&& tr);

  bool extension_respects_annotation
  (const std::vector<ZEvent>& trace, const ZAnnotation& annotation,
   const ZPartialOrder& mutated_po, const ZTrace& parent_trace) const;

  bool linearization_respects_annotation
  (const std::vector<ZEvent>& trace, const ZAnnotation& annotation,
   const ZPartialOrder& mutated_po, const ZTrace& parent_trace) const;

  bool global_variables_initialized_with_value_zero
  (const std::vector<ZEvent>& trace) const;


  /* *************************** */
  /* CONSTRUCTOR                 */
  /* *************************** */

 public:

  ~ZExplorer();

  ZExplorer(ZBuilderSC& tb);

 private:

  /* *************************** */
  /* STATISTICS                  */
  /* *************************** */

  // Number of fully executed traces
  unsigned executed_traces_full = 0;
  // Number of executed traces
  unsigned executed_traces = 1;
  // Number of times we used the interpreter to get a trace
  unsigned interpreter_used = 1;
  // Number of executed traces with some thread assume-blocked
  unsigned assume_blocked_thread = 0;
  // Number of 'full' traces ending in a deadlock
  unsigned executed_traces_full_deadlock = 0;
  // Number of annotated traces with no mutation
  // choices (eg all blocked by negative annotation)
  unsigned no_mut_choices = 0;
  // Number of mutations considered
  unsigned mutations_considered = 0;
  // Closure of mutated PO (already chrono-ordered) failed
  unsigned closure_failed = 0;
  // Closure of mutated PO (already chrono-ordered) succeeded
  unsigned closure_succeeded = 0;
  // Succeeded closure without any added edge
  unsigned closure_no_edge = 0;
  // Added edges during succeeded closures
  unsigned closure_edges = 0;
  // Outer-loop iterations of succeeded closures
  unsigned closure_iter = 0;
  // Total time spent on copying
  double time_copy = 0;
  // Total time spent on linearization
  double time_linearization = 0;
  // Total time spent on interpreting
  double time_interpreter = 0;
  // Total time spent on closure
  double time_closure = 0;
  // Total time spent on closure succ no edge
  double time_closure_no_edge = 0;
  // Linearization failed
  unsigned linearization_failed = 0;
  // Linearization succeeded
  unsigned linearization_succeeded = 0;
  // Linearization: num of parents and children, to estimate branching factor
  unsigned total_parents = 0;
  unsigned total_children = 0;
  // total linearization
  unsigned total_lin = 0;
  //avg branching factor
  double avg_branch = 0;
  // max branching factor
  double max_branch = 0;

};

#endif // __Z_EXPLORER_H__
