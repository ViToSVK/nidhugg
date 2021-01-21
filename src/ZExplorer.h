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

  std::unique_ptr<ZTrace> initial;

  bool info = false;

  class TraceExtension {
   public:
    TraceExtension() = default;
    TraceExtension
    (const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& ext_trace,
     const std::shared_ptr<ZGraph>& ext_graph,
     const std::shared_ptr<ZPartialOrder>& ext_po_full,
     const std::shared_ptr<ZPartialOrder>& ext_po_part,
     bool some_to_ann, bool assume_blocked);

    bool empty() const;

    std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>> trace;
    std::shared_ptr<ZGraph> graph;
    std::shared_ptr<ZPartialOrder> po_full;
    std::shared_ptr<ZPartialOrder> po_part;
    bool something_to_annotate;
    bool has_assume_blocked_thread;
    bool has_error;
  };

  /* *************************** */
  /* ALGORITHM                   */
  /* *************************** */

 public:

  mutable std::unordered_map<SymAddrSize,
                             std::map<ZEventID, ZTrace *>> parents;
  mutable std::unordered_map<
    SymAddrSize, std::map<ZEventID, std::set<ZTrace *>>> waitfor_negallowed;

  void process_backtrack_points
  (const ZPartialOrder& po_full, const ZEvent * write_lock) const;

  void process_backtrack_points_negallowed
  (const ZPartialOrder& po_full, const ZEvent * write_lock) const;

  bool explore();

  void print_stats() const;

 private:

  void add_backtrack_point
  (ZTrace * parent_trace, const CPid& cpid) const;

  bool try_add_backtrack_blocker
  (ZTrace * parent_trace, const CPid& cpid) const;

  bool try_add_backtrack_point
  (ZTrace * parent_trace, const CPid& cpid) const;

  bool explore_rec(ZTrace& ann_trace);

  bool mutate_read
  (ZTrace& ann_trace, const ZEvent *read, const ZAnn& mutation);

  const ZEvent * collect_lock_mutation
  (const ZTrace& ann_trace, const ZEvent *lock);

  bool mutate_lock
  (ZTrace& ann_trace, const ZEvent *lock, const ZEvent *unlock);

  bool close_po
  (ZTrace& ann_trace, const ZEvent *read_lock,
   ZAnnotation&& mutated_annotation, ZPartialOrder&& mutated_po,
   bool mutation_follows_current_trace);

  bool realize_mutation
  (ZTrace& parent_trace, const ZEvent *read_lock,
   ZAnnotation&& mutated_annotation, ZPartialOrder&& mutated_po,
   bool mutation_follows_current_trace);

  TraceExtension reuse_trace
  (const ZTrace& parent_trace, const ZEvent *read_lock,
   const ZAnnotation& mutated_annotation, ZPartialOrder&& mutated_po);

  TraceExtension extend_trace
  (std::vector<ZEvent>&& tr, ZPartialOrder&& mutated_po);

  bool early_stopping
  (const ZTrace& ann_trace,
   const std::map<const ZEvent *, std::set<ZAnn>>& read_mutations,
   const std::map<const ZEvent *, const ZEvent *>& lock_mutations);

  bool extension_respects_annotation
  (const std::vector<std::unique_ptr<ZEvent>>& trace, const ZAnnotation& annotation,
   const ZPartialOrder& mutated_po, const ZTrace& parent_trace) const;

  bool linearization_respects_annotation
  (const std::vector<ZEvent>& trace, const ZAnnotation& annotation,
   const ZPartialOrder& mutated_po, const ZTrace& parent_trace) const;

  bool global_variables_initialized_with_value_zero
  (const std::vector<std::unique_ptr<ZEvent>>& trace) const;


  /* *************************** */
  /* CONSTRUCTOR                 */
  /* *************************** */

 public:

  ZExplorer(ZBuilderSC& tb);

 private:

  /* *************************** */
  /* STATISTICS                  */
  /* *************************** */

  // Number of fully executed traces
  unsigned executed_traces_full = 0;
  // Number of fully+partially executed traces
  unsigned executed_traces = 1;
  // Number of times we used the interpreter to get a trace
  unsigned interpreter_used = 1;
  // Number of full/partial traces with some thread assume-blocked
  unsigned assume_blocked_thread = 0;
  // Number of 'full' traces ending in a deadlock
  unsigned executed_traces_full_deadlock = 0;
  // Early stopping considered and failed
  unsigned early_failed = 0;
  // Early stopping considered and succeeded
  unsigned early_succeeded = 0;
  // Number of reads with no mutation choices
  // (eg all blocked by negative annotation)
  unsigned no_mut_choices = 0;
  // Number of mutations considered
  unsigned mutations_considered = 0;
  // Closure of mutated POs failed
  unsigned closure_failed = 0;
  // Closure of mutated POs succeeded
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
  // Total time spent on early stopping
  double time_early = 0;
  // Linearization failed
  unsigned linearization_failed = 0;
  // Linearization succeeded
  unsigned linearization_succeeded = 0;
  // Linearization: num of parents and children, to estimate branching factor
  unsigned total_parents = 0;
  unsigned total_children = 0;
  // total linearization
  unsigned total_lin = 0;
  //avg2 branching factor
  double avg2_branch = 1.0;
  // max branching factor
  double max_branch = 1.0;
};

#endif // __Z_EXPLORER_H__
