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

#include <iostream>
#include <iomanip>

#include "VCExplorer.h"
#include "VCHelpers.h"
#include "VCDumps.cpp"

void VCExplorer::print_stats()
{
  std::setprecision(4);
  std::cout << "\n";
  std::cout << "Fully executed traces:            " << executed_traces_full << "\n";
  std::cout << "Fully+partially executed traces:  " << executed_traces << "\n";
  std::cout << "Interpreter used to get a trace:  " << interpreter_used << "\n";
  std::cout << "IntTr with assume-blocked thread: " << interpreter_assume_blocked_thread << "\n";
  std::cout << "Full traces ending in a deadlock: " << executed_traces_full_deadlock << "\n";
  std::cout << "Read-ordered partial orders:      " << read_ordered_pos << "\n";
  std::cout << "RoPOs with no mutation choices:   " << read_ordered_pos_no_mut_choices << "\n";
  std::cout << "Mutations considered:             " << mutations_considered << "\n";
  std::cout << "Closure after ordering failed:    " << cl_ordering_failed << "\n";
  std::cout << "Closure after ordering succeeded: " << cl_ordering_succeeded << "\n";
  std::cout << "Closure after mutation failed:    " << cl_mutation_failed << "\n";
  std::cout << "Closure after mutation succeeded: " << cl_mutation_succeeded << "\n";
  std::cout << "Total time spent on copying:      " << time_graphcopy << "\n";
  std::cout << "Total time spent on replaying:    " << time_replaying << "\n";
  std::cout << "Total time spent on mazurkiewicz: " << time_maz << "\n";
  std::cout << "Total time spent on closure:      " << time_closure << "\n";
  std::cout << "\n";

  // Change to false to test if assertions are on
  // To disable assertions (i.e. build as Release),
  // in src/Makefile add in CXXFLAGS this: -DNDEBUG
  assert(true && "RUN ON RELEASE");
}

/* *************************** */
/* EXPLORE                     */
/* *************************** */

bool VCExplorer::explore()
{
  while (!worklist.empty()) {
    // Get a VCTrace
    assert(!current.get());
    current = std::move(worklist.front());
    assert(!worklist.front().get());
    worklist.pop_front();
    mutationProducesMaxTrace.clear();

    //llvm::errs() << "********* TRACE *********\n";
    //current->annotation.dump();
    //current->graph.to_dot("");
    assert(traceRespectsAnnotation());

    // Get nodes available to be mutated
    auto nodesToMutate = current->graph.getNodesToMutate(current->annotation);
    assert(!nodesToMutate.empty());
    assert(current->graph.oneReadAndValueCausesMaxTrace == INT_MAX ||
           nodesToMutate.size() == 1);

    // Ordering of nodes to try mutations
    // Four VCDPOR version: mrl, mlr, rl, lr
    // 'm' means process of last mutation gets considered very first for new mutations
    // 'rl' means root process is considered before leaf processes
    // 'lr' considers leaf processes before root
    // leaves themselves are always ordered ascendingly by processID
    auto orderedNodesToMutate = orderNodesToMutate(nodesToMutate);

    std::vector<unsigned> processLengths = current->graph.getProcessLengths();

    // Get partial-order refinements that order extension events
    // Each refinement will be a candidate for possible mutations
    clock_t init = std::clock();
    std::list<PartialOrder> extendedPOs = orderingsAfterExtension();
    time_maz += (double)(clock() - init)/CLOCKS_PER_SEC;

    bool once = false;
    once_counted_full_trace = false;
    while (!extendedPOs.empty() || !once) {
      once = true;
      const PartialOrder& po =
        extendedPOs.empty() ? current->graph.getOriginal()
        : extendedPOs.front();
      bool haveOriginal = (po == current->graph.getOriginal());

      init = std::clock();
      auto withoutMutation = VCValClosure(current->graph, current->annotation);
      //llvm::errs() << "no-mutation-closure...";
      withoutMutation.valClose(po, nullptr, nullptr);
      //llvm::errs() << "done\n";
      time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

      if (!withoutMutation.closed) {
        // This ordering of extension events
        // makes the original annotation unrealizable
        // therefore no need to try any mutations
        ++cl_ordering_failed;
        assert(!haveOriginal && "current->trace is one witness");
        extendedPOs.front().first.reset();
        extendedPOs.front().second.reset();
        extendedPOs.pop_front();
        continue;
      }

      ++cl_ordering_succeeded;
      //llvm::errs() << "********* EXTENSION *********\n";
      //current->graph.to_dot(po, "");

      auto negativeWriteMazBranch = VCAnnotationNeg(current->negative);
      deadlockedExtension = true;

      // Try all possible nodes available to mutate
      for (auto ndit = orderedNodesToMutate.begin();
           ndit != orderedNodesToMutate.end(); ++ndit) {
        const Node * nd = *ndit;
        if (isRead(nd)) {
          deadlockedExtension = false;
          bool error = mutateRead(po, withoutMutation, negativeWriteMazBranch, nd);
          if (error) {
            assert(originalTB.error_trace);
            current.reset();
            worklist.clear();
            return error;
          }
          negativeWriteMazBranch.update(nd, processLengths);
        }
        else {
          assert(isLock(nd));
          bool error = mutateLock(po, withoutMutation, negativeWriteMazBranch, nd);
          if (error) {
            assert(originalTB.error_trace);
            current.reset();
            worklist.clear();
            return error;
          }
          negativeWriteMazBranch.update(nd, processLengths);
        }
      }

      // Done with this po
      if (deadlockedExtension) {
        // 'Full' trace ending in a deadlock
        ++executed_traces_full;
        ++executed_traces_full_deadlock;
      }
      if (!haveOriginal) {
        extendedPOs.front().first.reset();
        extendedPOs.front().second.reset();
        extendedPOs.pop_front();
      }
    } // end of loop for working with extension POs

    // Delete managed VCTrace
    current.reset();
  }

  return false;
}


/* *************************** */
/* ORDER NODES TO MUTATE       */
/* *************************** */

std::list<const Node *>
VCExplorer::orderNodesToMutate(std::unordered_set<const Node *>& nodesToMutate)
{
  auto result = std::list<const Node *>();
  auto& wrno = current->graph.scores_writeno;

  while (!nodesToMutate.empty()) {
    // Find (an arbitrary) maximal node (wrt our score)
    auto it = nodesToMutate.begin();
    assert(it != nodesToMutate.end());
    const Node *maxscore = *it;
    ++it;
    while (it != nodesToMutate.end()) {
      // Comparison
      bool yes = false;
      bool no = false;
      const Node *oth = *it;
      unsigned max_pid = maxscore->getProcessID();
      unsigned oth_pid = oth->getProcessID();
      // 1.) hides more writes than maxscore
      if (!yes && !no) {
        if (wrno.count(oth_pid) &&
            (!wrno.count(max_pid) || wrno[oth_pid] > wrno[max_pid]))
          yes = true;
        if (wrno.count(max_pid) &&
            (!wrno.count(oth_pid) || wrno[max_pid] > wrno[oth_pid]))
          no = true;
      }
      // 0.) process ID (root first)
      if (!yes && !no) {
        assert(max_pid != oth_pid);
        if (oth_pid == current->graph.starRoot())
          yes = true;
        else if (max_pid == current->graph.starRoot())
          no = true;
        else if (oth_pid < max_pid)
          yes = true;
        else
          no = true;
      }
      assert(yes != no);
      if (yes)
        maxscore = oth;
      ++it;
    }

    result.push_back(maxscore);
    nodesToMutate.erase(maxscore);
  }

  assert(nodesToMutate.empty());
  return result;
}


/* *************************** */
/* EXTENSION EVENTS ORDERINGS  */
/* *************************** */

std::list<PartialOrder> VCExplorer::orderingsAfterExtension()
{
  assert(current.get());
  current->graph.initWorklist();

  if (current->graph.lessThanTwoLeavesWithRorW())
    return current->graph.dumpDoneWorklist();

  // Go through all nonroot writes
  for (unsigned trace_idx = 0;
       trace_idx < current->trace.size(); ++trace_idx) {
    const VCEvent& ev = current->trace[trace_idx];
    if (current->graph.hasNodeWithEvent(ev)) {
      const Node *nd = current->graph.getNode(ev);
      if (isWrite(ev) && nd->getProcessID() != current->graph.starRoot()) {
        // Order with all nonroot annotated reads,
        // and with all nonroot writes such that:
        // 1) at least one is everGood
        // 2) both are observable
        current->graph.orderEventMaz(&ev, current->annotation, false,
                                     current->graph.getOriginal());
      }
    }
  }

  return current->graph.dumpDoneWorklist();
}

std::list<PartialOrder> VCExplorer::orderingsReadToBeMutated(const PartialOrder& po, const Node * nd)
{
  assert(isRead(nd));
  assert(current.get());
  assert(!current->annotation.defines(nd));
  current->graph.initWorklist();

  if (current->graph.lessThanTwoLeavesWithRorW())
    return current->graph.dumpDoneWorklist();

  // If the read is nonroot,
  // order it with all nonroot writes
  if (nd->getProcessID() != current->graph.starRoot())
    current->graph.orderEventMaz(nd->getEvent(), current->annotation,
                                 false, po);

  return current->graph.dumpDoneWorklist();
}

std::list<PartialOrder> VCExplorer::orderingsAfterMutationChoice
(const PartialOrder& po, const std::vector<const Node *> newEverGood)
{
  assert(current.get());
  // We have to create a po copy here since we try a mutation on it
  current->graph.initWorklist(po);

  if (current->graph.lessThanTwoLeavesWithRorW())
    return current->graph.dumpDoneWorklist();

  // After mutation choice, some nonrootwrites
  // become everGood, so they are from now required
  // to be ordered with other nonroot writes which
  // are notEverGood (if both are observable)
  for (auto& newEGnd : newEverGood)
    if (newEGnd->getProcessID() != current->graph.starRoot())
      current->graph.orderEventMaz(newEGnd->getEvent(), current->annotation, true, po);

  return current->graph.dumpDoneWorklist();
}

/* *************************** */
/* MUTATE READ                 */
/* *************************** */

bool VCExplorer::mutateRead(const PartialOrder& po, const VCValClosure& withoutMutation,
                            const VCAnnotationNeg& negativeWriteMazBranch, const Node *nd)
{
  assert(isRead(nd));
  // Orderings of the read just before mutating it
  clock_t init = std::clock();
  std::list<PartialOrder> readOrderedPOs = orderingsReadToBeMutated(po, nd);
  time_maz += (double)(clock() - init)/CLOCKS_PER_SEC;

  bool once = false;
  while (!readOrderedPOs.empty() || !once) {
    once = true;
    const PartialOrder& roPo =
      readOrderedPOs.empty() ? po : readOrderedPOs.front();
    bool haveArg = (roPo == po);
    ++read_ordered_pos;

    // We could do closure here to already rule out some po-s, but don't have to
    auto mutationCandidates =
      current->graph.getMutationCandidates(roPo, negativeWriteMazBranch, nd);
    if (mutationCandidates.empty()) {
      // All candidates are ruled out by negative annotation
      if (!haveArg) {
        readOrderedPOs.front().first.reset();
        readOrderedPOs.front().second.reset();
        readOrderedPOs.pop_front();
      }
      ++read_ordered_pos_no_mut_choices;
      continue;
    }

    int oneReadAndValueCausesMaxTrace = current->graph.oneReadAndValueCausesMaxTrace;
    for (auto& valpos_ann : mutationCandidates) {
      // llvm::errs() << valpos_ann.first.first << "_" << valpos_ann.first.second << "...";
      assert(valpos_ann.first.first != INT_MAX);
      if (oneReadAndValueCausesMaxTrace == valpos_ann.first.first) {
        // There is only one read to mutate, after it comes no read if it
        // sees this value, so we do not have to perform this mutation
        if (!once_counted_full_trace) {
          once_counted_full_trace = true;
          ++executed_traces_full;
        }
        continue;
      }

      if (mutationProducesMaxTrace.count(nd->getProcessID()) &&
          mutationProducesMaxTrace
          [nd->getProcessID()].count(valpos_ann.first.first)) {
        // This mutation was already done in a sibling Mazurkiewicz branch
        // And it produces a maximal trace with the same value function
        continue;
      }

      ++mutations_considered;
      // Collect writes that become newly everGood by performing this mutation
      // We will have to order them with conflicting notEverGood writes
      VCAnnotation mutatedAnnotation(current->annotation);
      auto newlyEverGoodVCIIDs = mutatedAnnotation.add(nd, valpos_ann.second);
      auto newlyEverGoodWrites = std::vector<const Node *>();
      for (auto& vciid : newlyEverGoodVCIIDs) {
        if (vciid.first != INT_MAX) {
          assert(isWrite(current->graph.getNode(vciid.first, vciid.second)));
          newlyEverGoodWrites.push_back(current->graph.getNode(vciid.first, vciid.second));
        }
      }
      // Sort the vector so that the execution is deterministic
      if (newlyEverGoodWrites.size() > 1) {
        auto comp = NodePtrComp();
        #ifndef NDEBUG
        for (auto it1 = newlyEverGoodWrites.begin();
             it1 != newlyEverGoodWrites.end(); ++it1)
          for (auto it2 = newlyEverGoodWrites.begin();
               it2 != it1; ++it2)
            assert(comp(*it1, *it2) || comp(*it2, *it1));
        #endif
        std::sort(newlyEverGoodWrites.begin(), newlyEverGoodWrites.end(), comp);
      }
      assert(mutatedAnnotation.size() == current->annotation.size() + 1);

      // Orderings of newly everGood writes
      clock_t init = std::clock();
      std::list<PartialOrder> afterMutationChoicePOs =
        orderingsAfterMutationChoice(roPo, newlyEverGoodWrites);
      assert(!afterMutationChoicePOs.empty());
      time_maz += (double)(clock() - init)/CLOCKS_PER_SEC;

      while (!afterMutationChoicePOs.empty()) {
        auto mutatedPo = std::move(afterMutationChoicePOs.front());
        assert(!afterMutationChoicePOs.front().first.get() &&
               !afterMutationChoicePOs.front().second.get());
        afterMutationChoicePOs.pop_front();

        // Closure after read+mutationChoice orderings and the actual mutation
        init = std::clock();
        auto withMutation = VCValClosure(withoutMutation);
        //llvm::errs() << "closure...";
        withMutation.valClose(mutatedPo, nd, &(valpos_ann.second));
        //llvm::errs() << "done\n";
        time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

        if (!withMutation.closed) {
          // The mutation on 'mutatedPo' failed
          ++cl_mutation_failed;
          // llvm::errs() << "FAILED\n";
          // current->graph.to_dot(mutatedPo, "");
          mutatedPo.first.reset();
          mutatedPo.second.reset();
          continue;
        }

        // The mutation on 'mutatedPo' succeeded
        ++cl_mutation_succeeded;
        // llvm::errs() << "SUCCEEDED\n";
        // current->graph.to_dot(mutatedPo, "");

        assert(mutatedAnnotation.size() == current->annotation.size() + 1);
        bool mutationFollowsCurrentTrace =
          (nd->getEvent()->value == valpos_ann.first.first);
        auto error_addedToWL = extendAndAdd(std::move(mutatedPo), mutatedAnnotation,
                                            negativeWriteMazBranch, nd->getProcessID(),
                                            mutationFollowsCurrentTrace);
        if (error_addedToWL.first)
          return true;
        if (!error_addedToWL.second) {
          // Annotating this read with this value on this trace produces a maximal trace
          // No need to repeat this mutation in a sibling Mazurkiewicz branch
          auto it = mutationProducesMaxTrace.find(nd->getProcessID());
          if (it == mutationProducesMaxTrace.end())
            mutationProducesMaxTrace.emplace_hint(it, nd->getProcessID(),
                                                  std::unordered_set<int>());
          assert(!mutationProducesMaxTrace
                 [nd->getProcessID()].count(valpos_ann.first.first));
          mutationProducesMaxTrace
            [nd->getProcessID()].insert(valpos_ann.first.first);
          // Clear the rest of POs with this mutation
          afterMutationChoicePOs.clear();
        }

      } // end of loop for newEverGoodOrdered partial order
    } // end of loop for mutation annotation
    if (!haveArg) {
        readOrderedPOs.front().first.reset();
        readOrderedPOs.front().second.reset();
        readOrderedPOs.pop_front();
    }
  } // end of loop for readOrdered partial order

  return false;
}

/* *************************** */
/* MUTATE LOCK                 */
/* *************************** */

bool VCExplorer::mutateLock(const PartialOrder& po, const VCValClosure& withoutMutation,
                            const VCAnnotationNeg& negativeWriteMazBranch, const Node *nd)
{
  assert(isLock(nd));
  auto lastLock = current->annotation.getLastLock(nd);

  if (!lastLock.first) {
    // This lock hasn't been touched before
    deadlockedExtension = false;
    if (negativeWriteMazBranch.forbidsInitialEvent(nd)) {
      // Negative annotation forbids initial unlock
      return false;
    }

    // Trivially realizable

    auto mutatedPo = PartialOrder(std::unique_ptr<ThreadPairsVclocks>
                                  (new ThreadPairsVclocks(*(po.first))),
                                  std::unique_ptr<ThreadPairsVclocks>
                                  (new ThreadPairsVclocks(*(po.second))));

    VCAnnotation mutatedAnnotation(current->annotation);
    mutatedAnnotation.setLastLock(nd);

    bool mutationFollowsCurrentTrace =
      (nd->getEvent()->observed_id == -1);
    auto error_addedToWL = extendAndAdd(std::move(mutatedPo), mutatedAnnotation,
                                        negativeWriteMazBranch, nd->getProcessID(),
                                        mutationFollowsCurrentTrace);
    return error_addedToWL.first;
  }

  // The lock has been touched before
  assert(lastLock.first);
  auto lastunlockit = current->graph.nodes_iterator(lastLock.second);
  #ifndef NDEBUG
  const Node * lastlocknd = *lastunlockit;
  assert(isLock(lastlocknd));
  #endif
  const Node * lastunlocknd = nullptr;

  while (!lastunlockit.atProcessEnd()) {
    ++lastunlockit;
    const VCEvent& cand = *((*lastunlockit)->getEvent());
    if (isUnlock(cand) && cand.ml == nd->getEvent()->ml) {
      assert(cand.ml == lastlocknd->getEvent()->ml);
      lastunlocknd = (*lastunlockit);
      break;
    }
  }

  if (!lastunlocknd) {
    // This lock is currently locked
    // Trivially unrealizable
    return false;
  }

  // This lock is currently unlocked by lastunlocknd
  assert(lastunlocknd);
  deadlockedExtension = false;

  if (negativeWriteMazBranch.forbids(nd, lastunlocknd)) {
    // Negative annotation forbids this unlock
    return false;
  }

  auto mutatedPo = PartialOrder(std::unique_ptr<ThreadPairsVclocks>
                                (new ThreadPairsVclocks(*(po.first))),
                                std::unique_ptr<ThreadPairsVclocks>
                                (new ThreadPairsVclocks(*(po.second))));

  clock_t init = std::clock();
  // Closure with lock-mutation
  auto withMutation = VCValClosure(withoutMutation);
  withMutation.valCloseLock(mutatedPo, nd, lastunlocknd);
  time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

  if (!withMutation.closed) {
    // The lock-mutation on 'mutatedPo' failed
    ++cl_mutation_failed;
    mutatedPo.first.reset();
    mutatedPo.second.reset();
    return false;
  }

  // The lock-mutation on 'mutatedPo' succeeded
  ++cl_mutation_succeeded;

  VCAnnotation mutatedAnnotation(current->annotation);
  mutatedAnnotation.setLastLock(nd);

  bool mutationFollowsCurrentTrace =
    (nd->getEvent()->observed_id == (int) lastunlocknd->getEvent()->id);
  auto error_addedToWL = extendAndAdd(std::move(mutatedPo), mutatedAnnotation,
                                      negativeWriteMazBranch, nd->getProcessID(),
                                      mutationFollowsCurrentTrace);
  return error_addedToWL.first;
}

/* *************************** */
/* EXTEND AND ADD              */
/* *************************** */

std::pair<bool, bool>
VCExplorer::extendAndAdd(PartialOrder&& mutatedPo,
                         const VCAnnotation& mutatedAnnotation,
                         const VCAnnotationNeg& negativeWriteMazBranch,
                         unsigned processMutationPreference,
                         bool mutationFollowsCurrentTrace)
{
  auto mutatedTrace = mutationFollowsCurrentTrace
    ?reuseTrace(mutatedAnnotation)
    :extendTrace(current->graph.linearize(mutatedPo, mutatedAnnotation));

  executed_traces++;
  if (mutatedTrace.hasError)
    return {true, false}; // Found an error
  if (mutatedTrace.hasAssumeBlockedThread) {
    // This recursion subtree of the algorithm will only
    // have traces that violate the same assume-condition
    interpreter_assume_blocked_thread++;
    //return false;
  }
  if (!mutatedTrace.somethingToAnnotate) {
    // Maximal trace, do not add to worklist
    //llvm::errs() << "********* FULL TRACE *********\n";
    //mutatedAnnotation.dump();
    //current->graph.to_dot(mutatedPo,"");
    executed_traces_full++;
    return {false, false};
  }

  clock_t init = std::clock();
  assert(mutatedPo.first.get() && mutatedPo.second.get());
  VCGraphVclock mutatedGraph(current->graph,       // base for graph
                             std::move(mutatedPo), // base for 'original' po
                             mutatedTrace.trace,   // to extend the graph
                             mutatedAnnotation);   // to extend the graph
  assert(!mutatedPo.first.get() && !mutatedPo.second.get());
  std::unique_ptr<VCTrace> mutatedVCTrace
    (new VCTrace(std::move(mutatedTrace.trace),
                 mutatedAnnotation,
                 negativeWriteMazBranch,
                 std::move(mutatedGraph),
                 processMutationPreference,
                 mutationFollowsCurrentTrace));
  assert(mutatedTrace.empty() && mutatedGraph.empty());
  time_graphcopy += (double)(clock() - init)/CLOCKS_PER_SEC;

  worklist.push_front(std::move(mutatedVCTrace));
  assert(!mutatedVCTrace.get());

  return {false, true};
}

/* *************************** */
/* REUSE TRACE                 */
/* *************************** */

VCExplorer::TraceExtension
VCExplorer::reuseTrace(const VCAnnotation& mutatedAnnotation)
{
  auto tr = std::vector<VCEvent>();
  tr.reserve(current->trace.size());
  bool somethingToAnnotate = false;


  for (const VCEvent& ev : current->trace) {
    if (!somethingToAnnotate) {
      const Node *nd = (current->graph.hasNodeWithEvent(ev))
        ?current->graph.getNode(ev):nullptr;
      if (isRead(ev) && (!nd || !mutatedAnnotation.defines(nd)))
        somethingToAnnotate = true;
      else if (isLock(ev)) {
        if (!nd ||
            (nd->getEventID() == current->graph[ nd->getProcessID() ].size() - 1 &&
             !mutatedAnnotation.isLastLock(nd)))
          somethingToAnnotate = true;
      }
    }

    tr.push_back(ev.copy(tr.size(), false));
  }

  return TraceExtension(std::move(tr), somethingToAnnotate);
}

/* *************************** */
/* EXTEND TRACE                */
/* *************************** */

VCExplorer::TraceExtension
VCExplorer::extendTrace(std::vector<VCEvent>&& tr)
{
  clock_t init = std::clock();
  VCTraceBuilder TB(originalTB.config, originalTB.M, std::move(tr));
  auto traceExtension = TraceExtension(TB.extendGivenTrace());
  time_replaying += (double)(clock() - init)/CLOCKS_PER_SEC;
  interpreter_used++;

  if (TB.has_error()) {
    // ERROR FOUND
    originalTB.error_trace = TB.get_trace();
    traceExtension.hasError = true;
  }

  if (TB.someThreadAssumeBlocked)
    traceExtension.hasAssumeBlockedThread = true;

  return traceExtension;
}

/* *************************** */
/* TRACE RESPECTS ANNOTATION   */
/* *************************** */

bool VCExplorer::traceRespectsAnnotation() const {
  std::unordered_map<SymAddrSize, VCIID> toObserve;
  std::unordered_set<SymAddrSize> dontObserve;
  for (unsigned i=0; i < current->trace.size(); ++i) {
    const VCEvent& ev = current->trace[i];
    if (isWrite(ev)) {
      if (!current->graph.hasNodeWithEvent(ev))
        dontObserve.insert(ev.ml);
      else
        toObserve[ev.ml] = VCIID(ev.pid, ev.event_order);
    }
    if (isRead(ev) && current->annotation.defines(ev.pid, ev.event_order)) {
      if (!current->reusedTrace && dontObserve.count(ev.ml)) {
        //current->graph.to_dot("");
        //current->annotation.dump();
        //current->graph.getNode(ev)->dump();
        //llvm::errs() << "OBSERVED WRITE NOT IN GRAPH\n";
        return false;
      }
      const auto& ann = current->annotation.getAnn(ev.pid, ev.event_order);
      if (ann.value != ev.value) {
        //current->graph.to_dot("");
        //current->annotation.dump();
        //current->graph.getNode(ev)->dump();
        //ann.dump();
        //llvm::errs() << "ANNVALUE: " << ann.value << "  EVENTVALUE: " << ev.value << "\n";
        //llvm::errs() << (current->reusedTrace?"REUSED":"FRESH") << " TRACE\n";
        return false;
      }
      if (!current->reusedTrace) {
        VCIID observed = (toObserve.count(ev.ml))?toObserve[ev.ml]:
          VCIID(INT_MAX, INT_MAX);
        bool observedIsGood = false;
        if (ann.goodLocal && (*(ann.goodLocal)) == observed)
          observedIsGood = true;
        else {
          for (const VCIID& gr : ann.goodRemote)
            if (gr == observed)
              observedIsGood = true;
        }
        if (!observedIsGood) {
          //current->graph.to_dot("");
          //current->annotation.dump();
          //current->graph.getNode(ev)->dump();
          //ann.dump();
          //llvm::errs() << "OBSERVED WRITE IS NOT GOOD: ["
          //             << observed.first << "][" << observed.second << "]\n";
          return false;
        }
      }
    }
  }

  return true;
}
