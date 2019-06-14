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

#ifndef __Z_EXPLORER_H__
#define __Z_EXPLORER_H__

#include <list>
#include <memory>
#include <time.h>

#include "ZBuilderTSO.h"
#include "ZTrace.h"
#include "ZClosure.h"


class ZExplorer {

  ZBuilderTSO& originalTB;

  const ZTrace *initial;

  bool info = true;

  class TraceExtension {
   public:
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

  bool exploreRec(const ZTrace& annTrace);

  bool mutateRead(const ZTrace& annTrace, const ZEvent *read);

  bool mutateLock(const ZTrace& annTrace, const ZEvent *lock);

  bool chronological
    (const ZTrace& annTrace, const ZEvent *readLock,
     ZAnnotation&& mutatedAnnotation, ZPartialOrder&& mutatedPO,
     bool mutationFollowsCurrentTrace);

  bool closePO
    (const ZTrace& annTrace, const ZEvent *readLock,
     ZAnnotation&& mutatedAnnotation, ZPartialOrder&& mutatedPO,
     bool mutationFollowsCurrentTrace);

  bool extendAndRecur
    (const ZTrace& parentTrace, const ZEvent *readLock,
     ZAnnotation&& mutatedAnnotation, ZPartialOrder&& mutatedPO,
     bool mutationFollowsCurrentTrace);

  TraceExtension reuseTrace(const ZTrace& parentTrace,
                            const ZAnnotation& mutatedAnnotation);

  TraceExtension extendTrace(std::vector<ZEvent>&& tr);

  bool respectsAnnotation(const std::vector<ZEvent>& trace,
                          const ZAnnotation& annotation,
                          const ZTrace& parentTrace) const;


  /* *************************** */
  /* CONSTRUCTOR                 */
  /* *************************** */

 public:

  ~ZExplorer();

  ZExplorer(ZBuilderTSO& tb);

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
  // Number of leaf-chrono-POs
  unsigned leaf_chrono_pos = 0;
  // Closure of mutated PO (already chrono-ordered) failed
  unsigned closure_failed = 0;
  // Closure of mutated PO (already chrono-ordered) succeeded
  unsigned closure_succeeded = 0;
  // Total time spent on copying
  double time_copy = 0;
  // Total time spent on linearization
  double time_linearization = 0;
  // Total time spent on interpreting
  double time_interpreter = 0;
  // Total time spent on chronological ordering
  double time_chrono = 0;
  // Total time spent on closure
  double time_closure = 0;

};

#endif // __Z_EXPLORER_H__
