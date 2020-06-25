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
#include "ZBuilderTSO.h"
#include "ZBuilderPSO.h"
#include "ZTrace.h"
#include "ZClosure.h"
#include "ZLinearization.h"


class ZExplorer {
 using BuffersT = std::unordered_map<std::vector<int>, std::unordered_map<SymAddrSize, std::list<ZEvent *>>>;
 public:
  TSOPSOTraceBuilder * original_tb = nullptr;
  const MemoryModel model;

 private:
  std::map<int, std::map<ZAnnotation, ZTrace>> schedules;
  std::map<int, std::set<ZAnnotation>> failed_schedules;
  std::map<int, std::set<ZAnnotation>> done_schedules;
  std::set<int> readlock_ids;

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
  // Linearization: num of parents and children, to estimate branching factor
  unsigned total_parents = 0;
  unsigned total_children = 0;
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
};

#endif // __Z_EXPLORER_H__
