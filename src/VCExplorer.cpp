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

void VCExplorer::print_stats()
{
  std::setprecision(4);
  std::cout << "\n";
  std::cout << "Fully executed traces:            " << executed_traces_full << "\n";
  std::cout << "Fully+partially executed traces:  " << executed_traces << "\n";
  std::cout << "F+P with assume-blocked thread:   " << executed_traces_assume_blocked_thread << "\n";
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

    //llvm::errs() << "********* TRACE *********\n";
    //current->annotation.dump();
    //current->graph.to_dot("");
    assert(traceRespectsAnnotation());

    // Get nodes available to be mutated
    auto nodesToMutate = current->graph.getNodesToMutate();
    for (auto it = nodesToMutate.begin(); it != nodesToMutate.end();) {
      const Node * nd = *it;
      assert(isRead(nd) || isLock(nd));
      if ((isRead(nd) && current->annotation.defines(nd)) ||
          (isLock(nd) && current->annotation.isLastLock(nd)))
        it = nodesToMutate.erase(it);
      else
        ++it;
    }

    #ifndef NDEBUG
    assert(!nodesToMutate.empty());
    assert(nodesToMutate.size() == current->unannot.size());
    for (auto& nd : nodesToMutate)
      assert(current->unannot.count(nd->getEvent()->iid.get_pid()));
    #endif

    // Ordering of nodes to try mutations
    // Four VCDPOR version: mrl, mlr, rl, lr
    // 'm' means process of last mutation gets considered very first for new mutations
    // 'rl' means root process is considered before leaf processes
    // 'lr' considers leaf processes before root
    // leaves themselves are always ordered ascendingly by processID
    auto orderedNodesToMutate = std::list<const Node *>();
    auto nonrootpid = std::vector<const Node *>(current->graph.size(), nullptr);
    const Node *pref = nullptr;
    const Node *root = nullptr;
    for (auto& ndtomut : nodesToMutate) {
      if (ndtomut->getProcessID() == current->processMutationPreference &&
          previous_mutation_process_first)
        pref = ndtomut;
      else if (ndtomut->getProcessID() == current->graph.starRoot())
        root = ndtomut;
      else
        nonrootpid.at(ndtomut->getProcessID()) = ndtomut;
    }

    // Nonroots are taken in the process-id-ascending fashion
    for (unsigned i=0; i<nonrootpid.size(); ++i)
      if (nonrootpid[i])
        orderedNodesToMutate.push_back(nonrootpid[i]);

    if (root) {
      if (root_before_nonroots)
        orderedNodesToMutate.push_front(root); // Root before nonroots
      else
        orderedNodesToMutate.push_back(root); // Root after nonroots
    }

    if (pref) {
      // Preference very first
      assert(previous_mutation_process_first);
      if (!root || pref != root)
        orderedNodesToMutate.push_front(pref);
    }

    std::vector<unsigned> processLengths = current->graph.getProcessLengths();

    // Get partial-order refinements that order extension events
    // Each refinement will be a candidate for possible mutations
    clock_t init = std::clock();
    std::list<PartialOrder> extendedPOs = orderingsAfterExtension();
    assert(!extendedPOs.empty() && "current->trace is one witness");
    time_maz += (double)(clock() - init)/CLOCKS_PER_SEC;

    while (!extendedPOs.empty()) {
      auto po = std::move(extendedPOs.front());
      assert(!extendedPOs.front().first.get() &&
             !extendedPOs.front().second.get());
      extendedPOs.pop_front();

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
        po.first.reset();
        po.second.reset();
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
      po.first.reset();
      po.second.reset();
    }

    // Delete managed VCTrace
    current.reset();
  }

  return false;
}

/* *************************** */
/* EXTENSION EVENTS ORDERINGS  */
/* *************************** */

std::list<PartialOrder> VCExplorer::orderingsAfterExtension()
{
  assert(current.get());
  current->graph.initWorklist();

  // Go through all nonroot writes
  for (unsigned trace_idx = 0;
       trace_idx < current->trace.size(); ++trace_idx) {
    const VCEvent& ev = current->trace[trace_idx];
    const Node *nd = current->graph.getNode(ev);
    if (isWrite(ev) && nd->getProcessID() != current->graph.starRoot()) {
      // Order with all nonroot annotated reads,
      // and with all nonroot writes such that:
      // 1) at least one is everGood
      // 2) both are observable
      current->graph.orderEventMaz(&ev, current->annotation, false);
    }
  }

  return current->graph.dumpDoneWorklist();
}

std::list<PartialOrder> VCExplorer::orderingsReadToBeMutated(const PartialOrder& po, const Node * nd)
{
  assert(isRead(nd));
  assert(current.get());
  assert(!current->annotation.defines(nd));
  current->graph.initWorklist(po);

  // If the read is nonroot,
  // order it with all nonroot writes
  if (nd->getProcessID() != current->graph.starRoot())
    current->graph.orderEventMaz(nd->getEvent(), current->annotation, false);

  return current->graph.dumpDoneWorklist();
}

std::list<PartialOrder> VCExplorer::orderingsAfterMutationChoice
(const PartialOrder& po, const std::vector<const Node *> newEverGood)
{
  assert(current.get());
  current->graph.initWorklist(po);

  // After mutation choice, some nonrootwrites
  // become everGood, so they are from now required
  // to be ordered with other nonroot writes which
  // are notEverGood (if both are observable)
  for (auto& newEGnd : newEverGood)
    if (newEGnd->getProcessID() != current->graph.starRoot())
      current->graph.orderEventMaz(newEGnd->getEvent(), current->annotation, true);

  return current->graph.dumpDoneWorklist();
}

/* *************************** */
/* MUTATE READ                 */
/* *************************** */

bool VCExplorer::mutateRead(const PartialOrder& po, const VCValClosure& withoutMutation,
                            const VCAnnotationNeg& negativeWriteMazBranch, const Node *nd)
{
  assert(isRead(nd));

  std::unordered_set<int> mutatedUnannot(current->unannot);
  assert(mutatedUnannot.count(nd->getEvent()->iid.get_pid()));
  mutatedUnannot.erase(nd->getEvent()->iid.get_pid());

  // Orderings of the read just before mutating it
  clock_t init = std::clock();
  std::list<PartialOrder> readOrderedPOs = orderingsReadToBeMutated(po, nd);
  assert(!readOrderedPOs.empty());
  time_maz += (double)(clock() - init)/CLOCKS_PER_SEC;

  while (!readOrderedPOs.empty()) {
    auto roPo = std::move(readOrderedPOs.front());
    assert(!readOrderedPOs.front().first.get() &&
           !readOrderedPOs.front().second.get());
    readOrderedPOs.pop_front();
    ++read_ordered_pos;

    // We could do closure here to already rule out some po-s, but don't have to
    auto mutationCandidates =
      current->graph.getMutationCandidates(roPo, negativeWriteMazBranch, nd);
    if (mutationCandidates.empty())
      ++read_ordered_pos_no_mut_choices;

    for (auto& valpos_ann : mutationCandidates) {

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
        auto error_addedToWL = extendAndAdd(std::move(mutatedPo), nullptr,
                                            mutatedUnannot, mutatedAnnotation,
                                            negativeWriteMazBranch, nd->getProcessID());
        if (error_addedToWL.first)
          return true;

      } // end of loop for newEverGoodOrdered partial order
    } // end of loop for mutation annotation
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

    std::unordered_set<int> mutatedUnannot(current->unannot);
    assert(mutatedUnannot.count(nd->getEvent()->iid.get_pid()));
    mutatedUnannot.erase(nd->getEvent()->iid.get_pid());

    auto mutatedPo = PartialOrder(std::unique_ptr<ThreadPairsVclocks>
                                  (new ThreadPairsVclocks(*(po.first))),
                                  std::unique_ptr<ThreadPairsVclocks>
                                  (new ThreadPairsVclocks(*(po.second))));

    VCAnnotation mutatedAnnotation(current->annotation);
    mutatedAnnotation.setLastLock(nd);

    auto mutatedLock = VCIID(nd->getProcessID(), nd->getEventID());
    auto error_addedToWL = extendAndAdd(std::move(mutatedPo), &mutatedLock,
                                        mutatedUnannot, mutatedAnnotation,
                                        negativeWriteMazBranch, nd->getProcessID());
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

  std::unordered_set<int> mutatedUnannot(current->unannot);
  assert(mutatedUnannot.count(nd->getEvent()->iid.get_pid()));
  mutatedUnannot.erase(nd->getEvent()->iid.get_pid());

  VCAnnotation mutatedAnnotation(current->annotation);
  mutatedAnnotation.setLastLock(nd);

  auto mutatedLock = VCIID(nd->getProcessID(), nd->getEventID());
  auto error_addedToWL = extendAndAdd(std::move(mutatedPo), &mutatedLock,
                                      mutatedUnannot, mutatedAnnotation,
                                      negativeWriteMazBranch, nd->getProcessID());
  return error_addedToWL.first;
}

/* *************************** */
/* EXTEND AND ADD              */
/* *************************** */

std::pair<bool, bool>
VCExplorer::extendAndAdd(PartialOrder&& mutatedPo,
                        const VCIID *mutatedLock,
                        const std::unordered_set<int>& mutatedUnannot,
                        const VCAnnotation& mutatedAnnotation,
                        const VCAnnotationNeg& negativeWriteMazBranch,
                        unsigned processMutationPreference)
{
  auto mutatedTrace =
    extendTrace(current->graph.linearize(mutatedPo, mutatedLock),
                mutatedUnannot);

  if (mutatedTrace.hasError)
    return {true, false}; // Found an error
  if (mutatedTrace.hasAssumeBlockedThread) {
    // This recursion subtree of the algorithm will only
    // have traces that violate the same assume-condition
    executed_traces_assume_blocked_thread++;
    //return false;
  }
  if (mutatedTrace.unannot.empty()) {
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
                             mutatedTrace.trace);  // to extend the graph
  assert(!mutatedPo.first.get() && !mutatedPo.second.get());
  std::unique_ptr<VCTrace> mutatedVCTrace
    (new VCTrace(std::move(mutatedTrace.trace),
                 mutatedAnnotation,
                 negativeWriteMazBranch,
                 std::move(mutatedGraph),
                 std::move(mutatedTrace.unannot),
                 processMutationPreference));
  assert(mutatedTrace.empty() && mutatedGraph.empty());
  time_graphcopy += (double)(clock() - init)/CLOCKS_PER_SEC;

  worklist.push_front(std::move(mutatedVCTrace));
  assert(!mutatedVCTrace.get());

  return {false, true};
}

/* *************************** */
/* EXTEND TRACE                */
/* *************************** */

VCExplorer::TraceExtension
VCExplorer::extendTrace(std::vector<VCEvent>&& tr,
                        const std::unordered_set<int>& unannot)
{
  clock_t init = std::clock();
  VCTraceBuilder TB(originalTB.config, originalTB.M, std::move(tr), unannot);
  auto traceExtension = TraceExtension(TB.extendGivenTrace());
  time_replaying += (double)(clock() - init)/CLOCKS_PER_SEC;
  executed_traces++;

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
  for (unsigned i=0; i < current->trace.size(); ++i) {
    const VCEvent& ev = current->trace[i];
    if (isRead(ev) && current->annotation.defines(ev.pid, ev.event_order)) {
      const auto& ann = current->annotation.getAnn(ev.pid, ev.event_order);
      if (ann.value != ev.value) {
        //current->graph.to_dot("");
        //current->annotation.dump();
        //current->graph.getNode(ev)->dump();
        //llvm::errs() << "ANNVALUE: " << ann.value << "  EVENTVALUE: " << ev.value << "\n";
        return false;
      }
      // TODO: check also if observes one of the good writes
    }
  }

  return true;
}
