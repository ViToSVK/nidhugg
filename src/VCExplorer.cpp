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

    // Get nodes available to be mutated
    auto nodesToMutate = current->graph.getNodesToMutate();

    assert(current->unannot.empty());
    for (auto& nd : nodesToMutate)
      current->unannot.insert(nd->getEvent()->iid.get_pid());

    for (auto it = nodesToMutate.begin(); it != nodesToMutate.end();) {
      const Node * nd = *it;
      assert(isRead(nd->getEvent()) || isLock(nd->getEvent()));
      if ((isRead(nd->getEvent()) && current->annotation.defines(nd)) ||
          (isLock(nd->getEvent()) && current->annotation.isLastLock(nd)))
        it = nodesToMutate.erase(it);
      else
        ++it;
    }

    if (nodesToMutate.empty()) {
      // Fully executed trace
      ++executed_traces_full;
      //llvm::errs() << "********* FULL TRACE *********\n";
      //current->annotation.dump();
      //current->graph.to_dot("");
      current.reset();
      continue;
    }

    if (!current->in_critical_section.empty()) {
      // Process in critical section, we mutate only on that process
      assert(current->in_critical_section.size() == 1);
      int cs_pid = current->in_critical_section.begin()->first;
      for (auto it = nodesToMutate.begin(); it != nodesToMutate.end(); ) {
        const Node * nd = *it;
        if (nd->getEvent()->iid.get_pid() == cs_pid)
          ++it;
        else
          it = nodesToMutate.erase(it);
      }
      assert(nodesToMutate.size() == 1);
    }

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
      if (!withoutMutation.closed) ++cl_ordering_failed;
      else ++cl_ordering_succeeded;

      if (!withoutMutation.closed) {
        // This ordering of extension events
        // makes the original annotation unrealizable
        // therefore no need to try any mutations
        po.first.reset();
        po.second.reset();
        continue;
      }

      //llvm::errs() << "********* EXTENSION *********\n";
      //current->graph.to_dot(po, "");

      auto negativeWriteMazBranch = VCAnnotationNeg(current->negative);

      // Try all possible nodes available to mutate
      for (auto ndit = orderedNodesToMutate.begin();
           ndit != orderedNodesToMutate.end(); ++ndit) {
        const Node * nd = *ndit;
        if (isRead(nd->getEvent())) {
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
          assert(isLock(nd->getEvent()));
          bool error = mutateLock(po, withoutMutation, negativeWriteMazBranch, nd);
          if (error) {
            assert(originalTB.error_trace);
            current.reset();
            worklist.clear();
            return error;
          }
        }
      }

      // Done with this po
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
  assert(isRead(nd->getEvent()));
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
(const PartialOrder& po, const std::unordered_set<const Node *> newEverGood)
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
  assert(isRead(nd->getEvent()));

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
      auto newlyEverGoodWrites = std::unordered_set<const Node *>();
      for (auto& vciid : newlyEverGoodVCIIDs) {
        if (vciid.first != INT_MAX) {
          assert(isWrite(current->graph.getNode(vciid.first, vciid.second)->getEvent()));
          newlyEverGoodWrites.insert(current->graph.getNode(vciid.first, vciid.second));
        }
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

        init = std::clock();
        auto mutatedTrace =
          extendTrace(current->graph.linearize(mutatedPo, nullptr),
                      mutatedUnannot);
        time_replaying += (double)(clock() - init)/CLOCKS_PER_SEC;
        if (mutatedTrace.first.empty())
          return true; // found an error
        assert(traceRespectsAnnotation(mutatedTrace.first, mutatedAnnotation));

        init = std::clock();
        VCGraphVclock mutatedGraph(current->graph,       // base for graph
                                   std::move(mutatedPo), // base for 'original' po
                                   mutatedTrace.first); // to extend the graph
        assert(!mutatedPo.first.get() && !mutatedPo.second.get());
        assert(mutatedAnnotation.size() == current->annotation.size() + 1);
        std::unique_ptr<VCTrace> mutatedVCTrace
          (new VCTrace(std::move(mutatedTrace.first),
                       mutatedAnnotation,
                       negativeWriteMazBranch,
                       std::move(mutatedGraph),
                       std::move(mutatedTrace.second),
                       nd->getProcessID()));
        assert(mutatedTrace.first.empty() &&
               mutatedGraph.empty() &&
               mutatedTrace.second.empty());
        time_graphcopy += (double)(clock() - init)/CLOCKS_PER_SEC;

        worklist.push_front(std::move(mutatedVCTrace));
        assert(!mutatedVCTrace.get());

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
  assert(isLock(nd->getEvent()));
  auto lastLock = current->annotation.getLastLock(nd);

  if (!lastLock.first) {
    // This lock hasn't been touched before
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
    clock_t init = std::clock();
    auto mutatedTrace =
      extendTrace(current->graph.linearize(mutatedPo, &mutatedLock),
                  mutatedUnannot);
    time_replaying += (double)(clock() - init)/CLOCKS_PER_SEC;
    if (mutatedTrace.first.empty())
      return true; // found an error
    assert(traceRespectsAnnotation(mutatedTrace.first, mutatedAnnotation));

    init = std::clock();
    VCGraphVclock mutatedGraph(current->graph,       // base for graph
                               std::move(mutatedPo), // base for 'original' po
                               mutatedTrace.first); // to extend the graph
    assert(!mutatedPo.first.get() && !mutatedPo.second.get());
    std::unique_ptr<VCTrace> mutatedVCTrace
      (new VCTrace(std::move(mutatedTrace.first),
                   mutatedAnnotation,
                   negativeWriteMazBranch,
                   std::move(mutatedGraph),
                   std::move(mutatedTrace.second),
                   nd->getProcessID()));
    assert(mutatedTrace.first.empty() &&
           mutatedGraph.empty() &&
           mutatedTrace.second.empty());
    time_graphcopy += (double)(clock() - init)/CLOCKS_PER_SEC;

    worklist.push_front(std::move(mutatedVCTrace));
    assert(!mutatedVCTrace.get());
    return false;
  }

  // The lock has been touched before
  assert(lastLock.first);
  auto lastunlockit = current->graph.nodes_iterator(lastLock.second);
  #ifndef NDEBUG
  const Node * lastlocknd = *lastunlockit;
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
  init = std::clock();
  auto mutatedTrace =
    extendTrace(current->graph.linearize(mutatedPo, &mutatedLock),
                mutatedUnannot);
  time_replaying += (double)(clock() - init)/CLOCKS_PER_SEC;
  if (mutatedTrace.first.empty())
    return true; // found an error
  assert(traceRespectsAnnotation(mutatedTrace.first, mutatedAnnotation));

  init = std::clock();
  VCGraphVclock mutatedGraph(current->graph,       // base for graph
                             std::move(mutatedPo), // base for 'original' po
                             mutatedTrace.first); // to extend the graph
  assert(!mutatedPo.first.get() && !mutatedPo.second.get());
  std::unique_ptr<VCTrace> mutatedVCTrace
    (new VCTrace(std::move(mutatedTrace.first),
                 mutatedAnnotation,
                 negativeWriteMazBranch,
                 std::move(mutatedGraph),
                 std::move(mutatedTrace.second),
                 nd->getProcessID()));
  assert(mutatedTrace.first.empty() &&
         mutatedGraph.empty() &&
         mutatedTrace.second.empty());
  time_graphcopy += (double)(clock() - init)/CLOCKS_PER_SEC;

  worklist.push_front(std::move(mutatedVCTrace));
  assert(!mutatedVCTrace.get());

  return false;
}

/* *************************** */
/* EXTEND TRACE                */
/* *************************** */

std::pair<std::vector<VCEvent>,
          std::unordered_map<int, int>>
VCExplorer::extendTrace(std::vector<VCEvent>&& tr,
                        const std::unordered_set<int>& unannot)
{
  VCTraceBuilder TB(originalTB.config, originalTB.M, std::move(tr), unannot);
  auto trace_cs = TB.extendGivenTrace();

  if (TB.has_error()) {
    // ERROR FOUND
    originalTB.error_trace = TB.get_trace();
    return {std::vector<VCEvent>(), std::unordered_map<int, int>()};
  }

  ++executed_traces;
  return trace_cs;
}

/* *************************** */
/* TRACE RESPECTS ANNOTATION   */
/* *************************** */

bool VCExplorer::traceRespectsAnnotation(const std::vector<VCEvent>& trace,
                                         const VCAnnotation& annotation) const {
  // TODO: check also if observes one of the good writes
  for (unsigned i=0; i < trace.size(); ++i) {
    const VCEvent& ev = trace[i];
    if (isRead(ev) && annotation.defines(ev.pid, ev.event_order)) {
      const auto& ann = annotation.getAnn(ev.pid, ev.event_order);
      if (ann.value != ev.value)
        return false;
      for (int j=i-1; j >= -1; --j) {
        if (j == -1) {

          if (ann.value != 0 || ev.value != 0)
            return false;
          if ((int) current->graph.starRoot() == ev.iid.get_pid() / 2 &&
              ann.loc != VCAnnotation::Loc::LOCAL)
            return false;
          if ((int) current->graph.starRoot() != ev.iid.get_pid() / 2 &&
              ann.loc != VCAnnotation::Loc::ANY)
            return false;
          break;

        } else {

          const VCEvent& wrev = trace[j];
          if (isWrite(wrev) && wrev.ml == ev.ml) {
            assert(wrev.value == ev.value);
            if (ann.value != wrev.value)
              return false;
            if ((int) current->graph.starRoot() == ev.iid.get_pid() / 2 &&
                ev.iid.get_pid() == wrev.iid.get_pid() &&
                ann.loc != VCAnnotation::Loc::LOCAL)
              return false;
            if ((int) current->graph.starRoot() == ev.iid.get_pid() / 2 &&
                ev.iid.get_pid() != wrev.iid.get_pid() &&
                ann.loc != VCAnnotation::Loc::REMOTE)
              return false;
            if ((int) current->graph.starRoot() != ev.iid.get_pid() / 2 &&
                ann.loc != VCAnnotation::Loc::ANY)
              return false;
            break;
          }

        }
      }
    }
  }

  return true;
}
