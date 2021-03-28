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

#include "ZBuilderSC.h"
#include "ZBuilderTSO.h"
#include "ZBuilderPSO.h"
#include "ZTrace.h"
#include "ZClosure.h"
#include "ZLinearization.h"
#include "ZLinNoclosure.h"
#include "ZLinNaive.h"


class ZExplorer {
 using BuffersT = std::unordered_map<std::vector<int>, std::unordered_map<SymAddrSize, std::list<ZEvent *>>>;
 public:
  TSOPSOTraceBuilder * original_tb = nullptr;
  const MemoryModel model;

 private:
  std::map<int, std::map<ZAnnotation, ZTrace>> schedules;
  std::map<int, std::set<ZAnnotation>> failed_schedules;
  std::map<int, std::set<ZAnnotation>> done_schedules;
  std::map<int, ZEventID> original_lock;
  std::set<int> readlock_ids;
  std::map<int, int> previous_lock_id;

  /* *************************** */
  /* CONSTRUCTOR                 */
  /* *************************** */

 public:
  ZExplorer(ZBuilderSC& tb);
  ZExplorer(ZBuilderTSO& tb);
  ZExplorer(ZBuilderPSO& tb);

  /* *************************** */
  /* STATISTICS                  */
  /* *************************** */

 private:
  // Number of fully executed traces
  unsigned executed_traces_full = 0;
  // Number of executed traces
  unsigned executed_traces = 0;
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
  // Closure of mutated PO failed
  unsigned closure_failed = 0;
  // Closure of mutated PO succeeded
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
  //avg2 branching factor
  double avg2_branch = 1.0;
  // max branching factor
  double max_branch = 1.0;
  //
  //
  //
  int lin_read_lower_bound = 10;
  int lin_perform_one_per = 5;
  int lin_goal = 5; //50;
  int lin_below_bound = 0;
  int lin_performed = 0;
  int lin_skipped = 0;
  int lin_not_even_rule1_succ = 0;
  std::vector<int> problematic_lin = std::vector<int>(1, -1);
  int repeat_runs = 10;
  //
  int max_allevents = 0;
  std::vector<int> no_allevents;
  std::vector<int> no_reads;
  std::vector<int> no_writes;
  std::vector<int> no_threads;
  std::vector<int> no_variables;
  std::vector<int> no_closure_rule23_edges;
  // times
  std::vector<double> t_rule1;
  std::vector<double> t_closure;
  std::vector<double> t_our_yescl_yesaux;
  std::vector<double> t_our_yescl_noaux;
  std::vector<double> t_our_nocl_yesaux;
  std::vector<double> t_our_nocl_noaux;
  std::vector<double> t_base_yescl_yesaux;
  std::vector<double> t_base_yescl_noaux;
  std::vector<double> t_base_nocl_yesaux;
  std::vector<double> t_base_nocl_noaux;
  // branching factors
  std::vector<double> br_our_yescl_yesaux;
  std::vector<double> br_our_yescl_noaux;
  std::vector<double> br_our_nocl_yesaux;
  std::vector<double> br_our_nocl_noaux;
  std::vector<double> br_base_yescl_yesaux;
  std::vector<double> br_base_yescl_noaux;
  std::vector<double> br_base_nocl_yesaux;
  std::vector<double> br_base_nocl_noaux;
  std::vector<int> par_our_yescl_yesaux;
  std::vector<int> par_our_yescl_noaux;
  std::vector<int> par_our_nocl_yesaux;
  std::vector<int> par_our_nocl_noaux;
  std::vector<int> par_base_yescl_yesaux;
  std::vector<int> par_base_yescl_noaux;
  std::vector<int> par_base_nocl_yesaux;
  std::vector<int> par_base_nocl_noaux;
  std::vector<int> ch_our_yescl_yesaux;
  std::vector<int> ch_our_yescl_noaux;
  std::vector<int> ch_our_nocl_yesaux;
  std::vector<int> ch_our_nocl_noaux;
  std::vector<int> ch_base_yescl_yesaux;
  std::vector<int> ch_base_yescl_noaux;
  std::vector<int> ch_base_nocl_yesaux;
  std::vector<int> ch_base_nocl_noaux;
  //
  //
  //
 public:
  void print_stats() const;

  /* *************************** */
  /* HELPERS                     */
  /* *************************** */

  void dump_schedules() const;
  void maintain_buffers(BuffersT& buffers, ZEvent * const ev, bool set_up_pointers) const;

  /* *************************** */
  /* ALGORITHM                   */
  /* *************************** */

 public:
  bool extend_and_explore(
    ZTrace& ann_trace, ZTraceExtension&& ext);

 private:
  bool explore(const ZTrace& ann_trace);

  void mutate(const ZTrace& ann_trace, const ZGraph& graph,
              const ZEvent * const readlock, const ZEventID& mutation);

  bool recur(const ZTrace& ann_trace);

  bool get_extension(ZTrace& ann_trace);

  bool linearization_respects_ann(
    const std::vector<ZEvent>& trace, const ZAnnotation& annotation,
    const ZGraph& graph, const ZTrace& parent_trace) const;

  void collect_linearization_stats(
    const ZGraph& graph);

  void linearization_experiments(
    const ZTrace& ann_trace,
    const ZAnnotation& annotation,
    const ZPartialOrder& closed_po,
    const ZPartialOrder& thread_order);
};

#endif // __Z_EXPLORER_H__
