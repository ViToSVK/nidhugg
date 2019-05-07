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

  std::list<std::unique_ptr<ZTrace>> worklist;

  std::unique_ptr<ZTrace> current;

  std::unordered_map<int, std::unordered_set<VCIID>>
    mutationProducesMaxTrace;

  class TraceExtension {
   public:
    TraceExtension(std::vector<ZEvent>&& extension, bool someToAnn)
      : trace(std::move(extension)),
      somethingToAnnotate(someToAnn),
      hasError(false), hasAssumeBlockedThread(false)
      {
        assert(extension.empty());
      }

    TraceExtension(std::pair<std::vector<ZEvent>&&, bool>&& extension_someToAnn)
      : TraceExtension(std::move(extension_someToAnn.first),
                       extension_someToAnn.second) {}

    bool empty() const {
      return (trace.empty());
    }

    std::vector<ZEvent> trace;
    bool somethingToAnnotate;
    bool hasError;
    bool hasAssumeBlockedThread;
  };

  /* *************************** */
  /* STATISTICS                  */
  /* *************************** */

  // Number of fully executed traces
  unsigned executed_traces_full = 0;
  // Number of executed traces
  unsigned executed_traces = 1;
  // Number of times we used the interpreter to get a trace
  unsigned interpreter_used = 1;
  // Number of interpreter-executed traces with some thread assume-blocked
  unsigned interpreter_assume_blocked_thread = 0;
  // Number of 'full' traces ending in a deadlock
  unsigned executed_traces_full_deadlock = 0;
  // Number of read-ordered partial orders
  unsigned read_ordered_pos = 0;
  // Number ofread-ordered partial orders with no mutation
  // choices (eg all blocked by negative annotation)
  unsigned read_ordered_pos_no_mut_choices = 0;
  // Number of mutations considered
  unsigned mutations_considered = 0;
  // Closure of extension-ordered po failed
  unsigned cl_ordering_failed = 0;
  // Closure of extension-ordered po succeeded
  unsigned cl_ordering_succeeded = 0;
  // Closure of mutation-chosen po failed
  unsigned cl_mutation_failed = 0;
  // Closure of mutation-chosen po succeeded
  unsigned cl_mutation_succeeded = 0;
  // Total time spent on copying
  double time_graphcopy = 0;
  // Total time spent on linearization and replaying
  double time_replaying = 0;
  // Total time spent on mazurkiewicz ordering
  double time_maz = 0;
  // Total time spent on closure
  double time_closure = 0;
  // Whether current extension-po ends in a deadlock
  bool deadlockedExtension;

  /* *************************** */
  /* ALGORITHM                   */
  /* *************************** */

  std::list<PartialOrder> orderingsAfterExtension();

  std::list<PartialOrder> orderingsReadToBeMutated(const PartialOrder& po, const Node * nd);

  std::list<PartialOrder> orderingsAfterMutationChoice
    (const PartialOrder& po, const std::vector<const Node *> newEverGood);

  bool mutateRead(const PartialOrder& po, const ZClosure& withoutMutation,
                  const ZAnnotationNeg& negativeWriteMazBranch, const Node *nd);

  bool mutateLock(const PartialOrder& po, const ZClosure& withoutMutation,
                  const ZAnnotationNeg& negativeWriteMazBranch, const Node *nd);

  std::pair<bool, bool> // <error?, added_into_worklist?>
    extendAndAdd(PartialOrder&& mutatedPo,
                 const ZAnnotation& mutatedAnnotation,
                 const ZAnnotationNeg& negativeWriteMazBranch,
                 unsigned processMutationPreference,
                 bool mutationFollowsCurrentTrace);

  TraceExtension reuseTrace(const ZAnnotation& mutatedAnnotation);

  TraceExtension extendTrace(std::vector<ZEvent>&& tr);

  bool traceRespectsAnnotation() const;

  bool previous_mutation_process_first = true;
  bool root_before_nonroots = true;

 public:

  bool explore();

  void print_stats();

  /* *************************** */
  /* CONSTRUCTOR                 */
  /* *************************** */

  ZExplorer(std::vector<ZEvent>&& initial_trace,
             bool somethingToAnnotate,
             ZBuilderTSO& tb,
             int star_root_index, bool p_m_p_f, bool r_b_n)
    : originalTB(tb),
    previous_mutation_process_first(p_m_p_f),
    root_before_nonroots(r_b_n)
    {
      if (!somethingToAnnotate)
        executed_traces_full = 1;
      else
        worklist.push_back(std::unique_ptr<ZTrace>(new ZTrace(std::move(initial_trace),
                                                                star_root_index)));
      if (tb.someThreadAssumeBlocked)
        interpreter_assume_blocked_thread = 1;
    }

};

#endif // __Z_EXPLORER_H__
