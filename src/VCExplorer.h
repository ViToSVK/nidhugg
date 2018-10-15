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

  /* *************************** */
  /* STATISTICS                  */
  /* *************************** */

  // Number of executed traces (1 because of the initial
  // trace that was obtained as a constructor argument)
  unsigned executed_traces = 1;
  // Number of fully executed traces
  unsigned executed_traces_full = 0;
  // Total time spent on closure
  double time_closure = 0;
  // Total time spent on copying
  double time_graphcopy = 0;
  // Total time spent on linearization and replaying
  double time_replaying = 0;
  // Total time spent on mazurkiewicz ordering
  double time_maz = 0;
  unsigned cl_ordering_failed = 0;
  unsigned cl_ordering_succeeded = 0;
  unsigned cl_mutation_failed = 0;
  unsigned cl_mutation_succeeded = 0;

  /* *************************** */
  /* ALGORITHM                   */
  /* *************************** */

  std::list<PartialOrder> extensionWritesOrderings();

  std::list<PartialOrder> readToBeMutatedOrderings(const PartialOrder& po, const Node * nd);

  std::list<PartialOrder> newlyObservableWritesOrderings
    (const PartialOrder& po, const std::unordered_set<const Node *> newobs);

  bool mutateRead(const PartialOrder& po, const VCValClosure& withoutMutation,
                  const VCAnnotationNeg& negativeWriteMazBranch, const Node *nd);

  bool mutateLock(const PartialOrder& po, const VCValClosure& withoutMutation, const Node *nd);

  std::pair<std::vector<VCEvent>,
    std::unordered_map<int, int>> extendTrace(std::vector<VCEvent>&& tr,
                                              const std::unordered_set<int>& unannot);

  bool traceRespectsAnnotation(const std::vector<VCEvent>& trace,
                               const VCAnnotation& annotation) const;

 public:

  bool explore();

  void print_stats();

  /* *************************** */
  /* CONSTRUCTOR                 */
  /* *************************** */

  VCExplorer(std::vector<VCEvent>&& trace, VCTraceBuilder& tb, int star_root_index)
    : originalTB(tb) {
    worklist.push_back(std::unique_ptr<VCTrace>(new VCTrace(std::move(trace), star_root_index)));
  }

};

#endif
