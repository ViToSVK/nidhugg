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

#include <iostream>
#include <iomanip>

#include "ZExplorer.h"
#include "ZHelpers.h"
#include "ZDumps.cpp"


ZExplorer::~ZExplorer()
{
  delete initial;
  initial = nullptr;
}


ZExplorer::ZExplorer
(std::vector<ZEvent>&& initial_trace,
 bool somethingToAnnotate,
 ZBuilderTSO& tb,
 int star_root_index)
  : originalTB(tb)
{
  if (tb.someThreadAssumeBlocked)
    interpreter_assume_blocked_thread = 1;

  if (!somethingToAnnotate)
    executed_traces_full = 1;
  else
    initial = new ZTrace(std::move(initial_trace), star_root_index);
}


void ZExplorer::print_stats()
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

bool ZExplorer::explore()
{
  if (initial)
    return exploreRec(*initial);
  return false;
}


bool ZExplorer::exploreRec(ZTrace& annTrace)
{
  assert(annTrace.respectsAnnotation());
  dumpTrace(annTrace.trace);
  annTrace.graph.dump();

  auto eventsToMutate = annTrace.getEventsToMutate();
  assert(!eventsToMutate.empty());
  for (auto& ev : eventsToMutate)
    ev->dump();

  // Unset this when any mutation succeeds
  annTrace.deadlocked = true;

  // Mutate the events
  for (const auto& ev : eventsToMutate) {
    if (isRead(ev)) {
      annTrace.deadlocked = false;
      bool error = mutateRead(annTrace, ev);
      if (error) {
        assert(originalTB.error_trace);
        return error;
      }
    }
    else {
      assert(isLock(ev));
      bool error = mutateLock(annTrace, ev);
      if (error) {
        assert(originalTB.error_trace);
        return error;
      }
    }
  }

  // Done with this recursion node
  if (annTrace.deadlocked) {
    // Not full trace ending in a deadlock
    // We count it as full
    ++executed_traces_full;
    ++executed_traces_full_deadlock;
  }

  return false;
}


/* *************************** */
/* MUTATE READ                 */
/* *************************** */

bool ZExplorer::mutateRead(const ZTrace& annTrace, const ZEvent *read)
{
  assert(isRead(read));

  auto obsCandidates = annTrace.getObsCandidates(read);
  if (obsCandidates.empty()) {
    // All candidates are ruled out by negative annotation
    ++read_ordered_pos_no_mut_choices;
    return false;
  }

  for (auto& observation : obsCandidates) {

    ++mutations_considered;
    ZAnnotation mutatedAnnotation(annTrace.annotation);
    mutatedAnnotation.add(read, observation);
    assert(mutatedAnnotation.size() == annTrace.annotation.size() + 1);
    auto mutatedPO = annTrace.graph.copyPO();

    auto init = std::clock();
    ZClosure preClosure(mutatedAnnotation, mutatedPO);
    preClosure.preClose(read, observation);
    time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

    mutatedAnnotation.dump();
    //bool mutationFollowsCurrentTrace =
    //  (nd->getEvent()->value == vciid_ann.second.value);
    // TODO call chronological
  }

  return false;
}


/* *************************** */
/* MUTATE LOCK                 */
/* *************************** */

bool ZExplorer::mutateLock(const ZTrace& annTrace, const ZEvent *lock)
{
  assert(isLock(lock));

  if (!annTrace.annotation.locationHasSomeLock(lock)) {
    // This lock hasn't been touched before
    annTrace.deadlocked = false;
    if (annTrace.negative.forbidsInitialEvent(lock)) {
      // Negative annotation forbids initial unlock
      return false;
    }

    // Trivially realizable
    ++mutations_considered;
    ZAnnotation mutatedAnnotation(annTrace.annotation);
    mutatedAnnotation.setLastLock(lock);
    auto mutatedPO = annTrace.graph.copyPO();

    mutatedAnnotation.dump();
    //bool mutationFollowsCurrentTrace =
    //(nd->getEvent()->observed_id == -1);
    // TODO call chronological
    return false;
  }

  // The lock has been touched before
  auto lastLockObs = annTrace.annotation.getLastLock(lock);
  auto lastUnlock = annTrace.graph.getBasis().getUnlockOfThisLock(lastLockObs);

  if (!lastUnlock) {
    // This lock is currently locked
    // Trivially unrealizable
    return false;
  }

  // This lock is currently unlocked by lastUnlock
  assert(lastUnlock && sameMl(lock, lastUnlock));
  annTrace.deadlocked = false;

  if (annTrace.negative.forbids(lock, lastUnlock)) {
    // Negative annotation forbids this unlock
    return false;
  }

  // Realizable
  ++mutations_considered;
  ZAnnotation mutatedAnnotation(annTrace.annotation);
  mutatedAnnotation.setLastLock(lock);
  auto mutatedPO = annTrace.graph.copyPO();

  auto init = std::clock();
  ZClosure preClosure(mutatedAnnotation, mutatedPO);
  preClosure.preClose(lock, lastUnlock);
  time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

  mutatedAnnotation.dump();
  //bool mutationFollowsCurrentTrace =
  //(nd->getEvent()->observed_id == (int) lastunlocknd->getEvent()->id);
  // TODO call chronological
  return false;
}



/*

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

  std::vector<unsigned> processLengths = current->graph.getProcessLengths();
  auto negativeWriteMazBranch = ZAnnotationNeg(current->negative);
  negativeWriteMazBranch.update(nd, processLengths);

  // Orderings of the read just before mutating it
  clock_t init = std::clock();
  std::list<PartialOrder> readOrderedPOs = orderingsReadToBeMutated(po, nd);
  time_maz += (double)(clock() - init)/CLOCKS_PER_SEC;


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
                                            std::unordered_set<VCIID>());
    assert(!mutationProducesMaxTrace
           [nd->getProcessID()].count(vciid_ann.first));
    mutationProducesMaxTrace
      [nd->getProcessID()].insert(vciid_ann.first);
    // Clear the rest of POs with this mutation
    afterMutationChoicePOs.clear();
  }

    auto error_addedToWL = extendAndAdd(std::move(mutatedPo), mutatedAnnotation,
                                        negativeWriteMazBranch, nd->getProcessID(),
                                        mutationFollowsCurrentTrace);
    return error_addedToWL.first;



/* *************************** */
/* EXTENSION EVENTS ORDERINGS  */
/* *************************** */
/*
std::list<PartialOrder> ZExplorer::orderingsAfterExtension()
{
  assert(current.get());
  current->graph.initWorklist();

  if (current->graph.lessThanTwoLeavesWithRorW())
    return current->graph.dumpDoneWorklist();

  // Go through all nonroot writes
  for (unsigned trace_idx = 0;
       trace_idx < current->trace.size(); ++trace_idx) {
    const ZEvent& ev = current->trace[trace_idx];
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

std::list<PartialOrder> ZExplorer::orderingsReadToBeMutated(const PartialOrder& po, const Node * nd)
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

std::list<PartialOrder> ZExplorer::orderingsAfterMutationChoice
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
/* EXTEND AND ADD              */
/* *************************** */
/*
std::pair<bool, bool>
ZExplorer::extendAndAdd(PartialOrder&& mutatedPo,
                         const ZAnnotation& mutatedAnnotation,
                         const ZAnnotationNeg& negativeWriteMazBranch,
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
  ZGraph mutatedGraph(current->graph,       // base for graph
                             std::move(mutatedPo), // base for 'original' po
                             mutatedTrace.trace,   // to extend the graph
                             mutatedAnnotation);   // to extend the graph
  assert(!mutatedPo.first.get() && !mutatedPo.second.get());
  std::unique_ptr<ZTrace> mutatedZTrace
    (new ZTrace(std::move(mutatedTrace.trace),
                 mutatedAnnotation,
                 negativeWriteMazBranch,
                 std::move(mutatedGraph),
                 processMutationPreference));
  assert(mutatedTrace.empty() && mutatedGraph.empty());
  time_graphcopy += (double)(clock() - init)/CLOCKS_PER_SEC;

  worklist.push_front(std::move(mutatedZTrace));
  assert(!mutatedZTrace.get());

  return {false, true};
}


/* *************************** */
/* REUSE TRACE                 */
/* *************************** */
/*
ZExplorer::TraceExtension
ZExplorer::reuseTrace(const ZAnnotation& mutatedAnnotation)
{
  auto tr = std::vector<ZEvent>();
  tr.reserve(current->trace.size());
  bool somethingToAnnotate = false;


  for (const ZEvent& ev : current->trace) {
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
/*
ZExplorer::TraceExtension
ZExplorer::extendTrace(std::vector<ZEvent>&& tr)
{
  clock_t init = std::clock();
  ZBuilderTSO TB(originalTB.config, originalTB.M, std::move(tr));
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
/*
bool ZExplorer::traceRespectsAnnotation() const {
  for (unsigned i=0; i < current->trace.size(); ++i) {
    const ZEvent& ev = current->trace[i];
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
*/
