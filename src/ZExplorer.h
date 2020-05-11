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

  TSOPSOTraceBuilder * originalTB = nullptr;

  ZTrace * initial;

  bool info = false;
  const bool tso;
  // SC flag is passed only to explorer, so it knows which builder to invoke
  // during extensions; ZGraph is not passed this flag, it handles SC as TSO
  const bool sc_flag;

  class TraceExtension {
   public:
    TraceExtension() = default;
    TraceExtension(std::vector<ZEvent>&& extension,
                 bool someToAnn, bool assumeBlocked);

    TraceExtension(std::pair<std::vector<ZEvent>&&, bool>&& extension_someToAnn);

    bool empty() const { return (trace.empty()); }

    std::vector<ZEvent> trace;
    bool somethingToAnnotate;
    bool hasAssumeBlockedThread;
    bool hasError;
  };


  /* *************************** */
  /* ALGORITHM                   */
  /* *************************** */

 public:

  bool explore();

  void print_stats() const;

 private:

  bool exploreRec(ZTrace& annTrace);

  bool mutateRead(const ZTrace& annTrace, const ZEvent *read);

  bool mutateLock(const ZTrace& annTrace, const ZEvent *lock);

  bool closePO
    (const ZTrace& annTrace, const ZEvent *readLock,
     ZAnnotation&& mutatedAnnotation, ZPartialOrder&& mutatedPO);

  bool extendAndRecur
    (const ZTrace& parentTrace, ZAnnotation&& mutatedAnnotation,
     ZPartialOrder&& mutatedPO);

  TraceExtension extendTrace(std::vector<ZEvent>&& tr);

  bool respectsAnnotation(const std::vector<ZEvent>& trace,
                          const ZAnnotation& annotation,
                          const ZPartialOrder& mutatedPO,
                          const ZTrace& parentTrace) const;

  bool linearizationRespectsAnn(const std::vector<ZEvent>& trace,
                                const ZAnnotation& annotation,
                                const ZPartialOrder& mutatedPO,
                                const ZTrace& parentTrace) const;


  /* *************************** */
  /* CONSTRUCTOR                 */
  /* *************************** */

 public:

  ~ZExplorer();

  ZExplorer(ZBuilderSC& tb);
  ZExplorer(ZBuilderTSO& tb);
  ZExplorer(ZBuilderPSO& tb);

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
  // Linearization: num of parents and children, to estimate branching factor
  unsigned total_parents = 0;
  unsigned total_children = 0;

};

#endif // __Z_EXPLORER_H__
