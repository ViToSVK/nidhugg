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

#ifndef __VC_EXPLORER_H__
#define __VC_EXPLORER_H__

#include <list>
#include <memory>
#include <time.h>

#include "VCTraceBuilder.h"
#include "VCTrace.h"
#include "VCValClosure.h"

class VCExplorer {

  VCTraceBuilder& originalTB;

  std::list<std::unique_ptr<VCTrace>> worklist;

  std::unique_ptr<VCTrace> current;

  class TraceExtension {
   public:
    TraceExtension(std::pair<std::vector<VCEvent>,
                             std::unordered_set<int>>&& extension_unannot)
      : trace(std::move(extension_unannot.first)),
      unannot(std::move(extension_unannot.second)),
      hasError(false), hasAssumeBlockedThread(false)
      {
        assert(extension_unannot.first.empty());
        assert(extension_unannot.second.empty());
      }

    bool empty() const {
      return (trace.empty() && unannot.empty());
    }

    std::vector<VCEvent> trace;
    std::unordered_set<int> unannot;
    bool hasError;
    bool hasAssumeBlockedThread;
  };

  /* *************************** */
  /* STATISTICS                  */
  /* *************************** */

  // Number of executed traces (1 because of the initial
  // trace that was obtained as a constructor argument)
  // This is the number of times we used the interpreter
  unsigned executed_traces = 1;
  // Number of fully executed traces
  unsigned executed_traces_full = 0;
  // Number of 'full' traces ending in a deadlock
  unsigned executed_traces_full_deadlock = 0;
  // Number of executed traces with some thread assume-blocked
  unsigned executed_traces_assume_blocked_thread = 0;
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
  // Total time spent on closure
  double time_closure = 0;
  // Total time spent on copying
  double time_graphcopy = 0;
  // Total time spent on linearization and replaying
  double time_replaying = 0;
  // Total time spent on mazurkiewicz ordering
  double time_maz = 0;
  // Whether current extension-po ends in a deadlock
  bool deadlockedExtension;

  /* *************************** */
  /* ALGORITHM                   */
  /* *************************** */

  std::list<PartialOrder> orderingsAfterExtension();

  std::list<PartialOrder> orderingsReadToBeMutated(const PartialOrder& po, const Node * nd);

  std::list<PartialOrder> orderingsAfterMutationChoice
    (const PartialOrder& po, const std::vector<const Node *> newEverGood);

  bool mutateRead(const PartialOrder& po, const VCValClosure& withoutMutation,
                  const VCAnnotationNeg& negativeWriteMazBranch, const Node *nd);

  bool mutateLock(const PartialOrder& po, const VCValClosure& withoutMutation,
                  const VCAnnotationNeg& negativeWriteMazBranch, const Node *nd);

  std::pair<bool, bool> // <error?, added_into_worklist?>
    extendAndAdd(PartialOrder&& mutatedPo,
                 const VCIID *mutatedLock,
                 const std::unordered_set<int>& mutatedUnannot,
                 const VCAnnotation& mutatedAnnotation,
                 const VCAnnotationNeg& negativeWriteMazBranch,
                 unsigned processMutationPreference);

  TraceExtension extendTrace(std::vector<VCEvent>&& tr,
                             const std::unordered_set<int>& unannot);

  bool traceRespectsAnnotation() const;

  bool previous_mutation_process_first = true;
  bool root_before_nonroots = true;

 public:

  bool explore();

  void print_stats();

  /* *************************** */
  /* CONSTRUCTOR                 */
  /* *************************** */

  VCExplorer(std::vector<VCEvent>&& initial_trace,
             std::unordered_set<int>&& initial_unannot,
             VCTraceBuilder& tb,
             int star_root_index, bool p_m_p_f, bool r_b_n)
    : originalTB(tb),
    previous_mutation_process_first(p_m_p_f),
    root_before_nonroots(r_b_n)
    {
      if (!initial_unannot.empty())
        worklist.push_back(std::unique_ptr<VCTrace>(new VCTrace(std::move(initial_trace),
                                                                std::move(initial_unannot),
                                                                star_root_index)));
      if (tb.someThreadAssumeBlocked)
        executed_traces_assume_blocked_thread = 1;
    }

};

#endif
