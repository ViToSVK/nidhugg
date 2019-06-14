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


ZExplorer::TraceExtension::TraceExtension
(std::vector<ZEvent>&& extension, bool someToAnn, bool assumeBlocked)
  : trace(std::move(extension)), somethingToAnnotate(someToAnn),
    hasAssumeBlockedThread(assumeBlocked), hasError(false)
{
  assert(extension.empty());
}


ZExplorer::TraceExtension::TraceExtension
(std::pair<std::vector<ZEvent>&&, bool>&& extension_someToAnn)
  : TraceExtension(std::move(extension_someToAnn.first),
                   extension_someToAnn.second, false) {}


ZExplorer::~ZExplorer()
{
  delete initial;
  initial = nullptr;
}


ZExplorer::ZExplorer(ZBuilderTSO& tb)
  : originalTB(tb)
{
  if (tb.someThreadAssumeBlocked)
    assume_blocked_thread = 1;

  if (tb.somethingToAnnotate.empty())
    executed_traces_full = 1;
  else
    initial = new ZTrace(std::move(tb.prefix), tb.star_root_index,
                         tb.someThreadAssumeBlocked);
}


void ZExplorer::print_stats() const
{
  std::setprecision(4);
  std::cout << "\n";
  std::cout << "Fully executed traces:             " << executed_traces_full << "\n";
  std::cout << "Fully+partially executed traces:   " << executed_traces << "\n";
  std::cout << "Interpreter used to get a trace:   " << interpreter_used << "\n";
  std::cout << "Traces with assume-blocked thread: " << assume_blocked_thread << "\n";
  std::cout << "Full traces ending in a deadlock:  " << executed_traces_full_deadlock << "\n";
  std::cout << "Traces with no mutation choices:   " << no_mut_choices << "\n";
  std::cout << "Mutations considered:              " << mutations_considered << "\n";
  std::cout << "Leaf-chronological-POs:            " << leaf_chrono_pos << "\n";
  std::cout << "Closure failed:                    " << closure_failed << "\n";
  std::cout << "Closure succeeded:                 " << closure_succeeded << "\n";
  std::cout << "Time spent on copying:             " << time_copy << "\n";
  std::cout << "Time spent on linearization:       " << time_linearization << "\n";
  std::cout << "Time spent on interpreting:        " << time_interpreter << "\n";
  std::cout << "Time spent on chronological:       " << time_chrono << "\n";
  std::cout << "Time spent on closure:             " << time_closure << "\n";
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


bool ZExplorer::exploreRec(const ZTrace& annTrace)
{
  if (info) {
    //dumpTrace(annTrace.trace);
    annTrace.graph.dump();
    llvm::errs() << "-------------------------\n\n";
  }

  auto eventsToMutate = annTrace.getEventsToMutate();
  assert(!eventsToMutate.empty());

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
    ++no_mut_choices;
    return false;
  }

  for (auto& observation : obsCandidates) {
    ++mutations_considered;

    ZAnnotation mutatedAnnotation(annTrace.annotation);
    mutatedAnnotation.add(read, observation);
    assert(mutatedAnnotation.size() == annTrace.annotation.size() + 1);
    auto mutatedPO = annTrace.graph.copyPO();

    if (info) {
      llvm::errs() << "Read to mutate:\n";
      read->dump();
      llvm::errs() << "Observation:\n" << observation.to_string() << "\n";
      mutatedAnnotation.dump();
    }

    auto init = std::clock();
    ZClosure preClosure(mutatedAnnotation, mutatedPO);
    preClosure.preClose(read, observation);
    time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

    bool mutationFollowsCurrentTrace = false;
    if (observation.thr == INT_MAX) {
      mutationFollowsCurrentTrace = //VC (read->value == 0);
        (read->observed_trace_id == -1);
    } else {
      const ZEvent *obsB = annTrace.getEvent(observation);
      assert(isWriteB(obsB) && isWriteM(obsB->write_other_ptr) &&
             sameMl(read, obsB) && sameMl(read, obsB->write_other_ptr));
      mutationFollowsCurrentTrace = //VC (read->value == obsB->value);
        (read->observed_trace_id == obsB->traceID() ||
         read->observed_trace_id == obsB->write_other_ptr->traceID());
    }

    bool error = chronological
      (annTrace, read, std::move(mutatedAnnotation),
       std::move(mutatedPO), mutationFollowsCurrentTrace);

    if (error) {
      assert(originalTB.error_trace);
      return error;
    }
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

    if (info) {
      llvm::errs() << "Lock to mutate:\n";
      lock->dump();
      llvm::errs() << "Not locked before\n";
      mutatedAnnotation.dump();
    }

    bool mutationFollowsCurrentTrace =
      (lock->observed_trace_id == -1);

    return chronological
      (annTrace, lock, std::move(mutatedAnnotation),
       std::move(mutatedPO), mutationFollowsCurrentTrace);
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
  assert(lastUnlock && isUnlock(lastUnlock) &&
         sameMl(lock, lastUnlock));
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

  if (info) {
    llvm::errs() << "Lock to mutate:\n";
    lock->dump();
    llvm::errs() << "Currently unlocked by:\n";
    lastUnlock->dump();
    mutatedAnnotation.dump();
  }

  auto init = std::clock();
  ZClosure preClosure(mutatedAnnotation, mutatedPO);
  preClosure.preClose(lock, lastUnlock);
  time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

  bool mutationFollowsCurrentTrace =
    (lock->observed_trace_id == lastUnlock->traceID());

  return chronological
    (annTrace, lock, std::move(mutatedAnnotation),
     std::move(mutatedPO), mutationFollowsCurrentTrace);
}


/* *************************** */
/* CHRONOLOGICAL               */
/* *************************** */

bool ZExplorer::chronological
(const ZTrace& annTrace, const ZEvent *readLock,
 ZAnnotation&& mutatedAnnotation, ZPartialOrder&& mutatedPO,
 bool mutationFollowsCurrentTrace)
{
  assert(isRead(readLock) || isLock(readLock));
  int chronoThreads = annTrace.graph.getCache().chrono.size();
  if (isRead(readLock) && !annTrace.isRoot(readLock) &&
      !annTrace.graph.getCache().chrono.count(readLock->threadID()))
    ++chronoThreads;

  if (chronoThreads < 2) {
    // No chronological orderings are needed,
    // we proceed with just one mutatedPO
    ++leaf_chrono_pos;

    // TODO Optimization:
    // If we want r to observe remote goodwm
    // and in trace r is before goodwm
    // we can just delay thread of r
    // such that in new trace
    // badlocalwm -> goodwm -> r

    return closePO
      (annTrace, readLock, std::move(mutatedAnnotation),
       std::move(mutatedPO), mutationFollowsCurrentTrace);
  }

  auto toOrder = annTrace.graph.chronoOrderPairs
    ((isRead(readLock) && !annTrace.isRoot(readLock))
     ? readLock : nullptr);

  if (toOrder.empty()) {
    // Nothing to order
    ++leaf_chrono_pos;
    return closePO
      (annTrace, readLock, std::move(mutatedAnnotation),
       std::move(mutatedPO), mutationFollowsCurrentTrace);
  }

  // First - partial order
  // Second - check toOrder pairs starting from this one
  std::list<std::pair<ZPartialOrder, unsigned>> worklist;
  worklist.emplace_front(std::move(mutatedPO), 0);
  assert(mutatedPO.empty() && !worklist.empty() &&
         !worklist.front().first.empty());

  while (!worklist.empty()) {
    // TODO Optimization: check if this mutation leads to
    // full trace, if yes you can break from this loop

    // Recursively process one chronoPO branch
    std::pair<ZPartialOrder, unsigned> current
      (std::move(worklist.front()));
    assert(!current.first.empty() && !worklist.empty() &&
           worklist.front().first.empty());
    worklist.pop_front();

    while (current.second < toOrder.size()) {
      auto ev1 = toOrder[current.second].first;
      auto ev2 = toOrder[current.second].second;
      // Order pair ev1-ev2
      // TODO: add *one-read*-or-*both-mws-observable-in-po* condition
      if (!current.first.areOrdered(ev1, ev2)) {
        // Create two cases with these orderings
        worklist.emplace_front(ZPartialOrder(current.first), // copy
                               current.second + 1);
        // Handle current: ev1 -> ev2
        current.first.addEdge(ev1, ev2);
        // Handle otherorder: ev2 -> ev1
        worklist.front().first.addEdge(ev2, ev1);
      }
      current.second++;
    }
    assert(current.second == toOrder.size());
    #ifndef NDEBUG
    for (auto& pair : toOrder)
      current.first.areOrdered(pair.first, pair.second);
      // or-*one-read*-or-*both-mws-observable-in-po* condition
    #endif
    ++leaf_chrono_pos;

    bool error = closePO
      (annTrace, readLock,
       worklist.empty()
       ? std::move(mutatedAnnotation) // move
       : ZAnnotation(mutatedAnnotation), // copy
       std::move(current.first), mutationFollowsCurrentTrace);

    if (error) {
      assert(originalTB.error_trace);
      return error;
    }
  }
  worklist.clear();

  return false;
}


/* *************************** */
/* CLOSE PO                    */
/* *************************** */

bool ZExplorer::closePO
(const ZTrace& annTrace, const ZEvent *readLock,
 ZAnnotation&& mutatedAnnotation, ZPartialOrder&& mutatedPO,
 bool mutationFollowsCurrentTrace)
{
  auto init = std::clock();
  ZClosure closure(mutatedAnnotation, mutatedPO);
  bool closed = closure.close
    (isLock(readLock) ? nullptr : readLock);
  time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

  if (!closed) {
    if (info) llvm::errs() << "Closure failed\n\n-----\n\n";
    ++closure_failed;
    return false;
  }

  if (info) llvm::errs() << "Closure succeeded\n\n-----\n\n";
  ++closure_succeeded;

  return extendAndRecur
    (annTrace, readLock, std::move(mutatedAnnotation),
     std::move(mutatedPO), mutationFollowsCurrentTrace);
}


/* *************************** */
/* EXTEND AND RECUR            */
/* *************************** */

bool ZExplorer::extendAndRecur
(const ZTrace& parentTrace, const ZEvent *readLock,
 ZAnnotation&& mutatedAnnotation, ZPartialOrder&& mutatedPO,
 bool mutationFollowsCurrentTrace)
{
  assert(isRead(readLock) || isLock(readLock));
  auto mutatedTrace = mutationFollowsCurrentTrace
    ? reuseTrace(parentTrace, mutatedAnnotation)
    : extendTrace(parentTrace.graph.linearize(mutatedPO, mutatedAnnotation));

  assert(respectsAnnotation(mutatedTrace.trace, mutatedAnnotation, parentTrace));

  executed_traces++;
  if (mutatedTrace.hasError) {
    assert(!mutationFollowsCurrentTrace);
    return true; // Found an error
  }
  if (mutatedTrace.hasAssumeBlockedThread) {
    // This recursion subtree of the algorithm will only
    // have traces that violate the same assume-condition
    assume_blocked_thread++;
    //return false;
  }
  if (!mutatedTrace.somethingToAnnotate) {
    // Maximal trace
    executed_traces_full++;
    // TODO Optimization:
    // Note in explorer that mutation produces max-trace
    return false;
  }

  clock_t init = std::clock();
  ZAnnotationNeg mutatedNeg(parentTrace.negative);
  if (isRead(readLock)) {
    mutatedNeg.update
      (readLock, parentTrace.graph.getBasis().real_sizes());
  }
  assert(!mutatedTrace.empty() && !mutatedPO.empty() &&
         !mutatedAnnotation.empty() && !mutatedNeg.empty());
  ZTrace mutatedZTrace
    (parentTrace,
     std::move(mutatedTrace.trace),
     std::move(mutatedAnnotation),
     std::move(mutatedNeg),
     std::move(mutatedPO),
     mutatedTrace.hasAssumeBlockedThread);
  assert(mutatedTrace.empty() && mutatedPO.empty() &&
         mutatedAnnotation.empty() && mutatedNeg.empty());
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;

  return exploreRec(mutatedZTrace);
}


/* *************************** */
/* REUSE TRACE                 */
/* *************************** */

ZExplorer::TraceExtension
ZExplorer::reuseTrace
(const ZTrace& parentTrace, const ZAnnotation& mutatedAnnotation)
{
  std::vector<ZEvent> tr;
  tr.reserve(parentTrace.trace.size());
  bool somethingToAnnotate = false;

  for (const ZEvent& ev : parentTrace.trace) {
    if (!somethingToAnnotate) {
      if (isRead(ev) && (!mutatedAnnotation.defines(&ev)))
        somethingToAnnotate = true;
      else if (isLock(ev)) {
        if (!parentTrace.graph.getBasis().hasEvent(&ev))
          somethingToAnnotate = true;
        else if (ev.eventID() == // last in its thraux
                 (parentTrace.graph.getBasis())(ev.threadID(), ev.auxID()).size() - 1
                 && !mutatedAnnotation.isLastLock(&ev))
          somethingToAnnotate = true;
      }
    }

    tr.push_back(ev.copy(tr.size(), false));
  }

  return TraceExtension(std::move(tr), somethingToAnnotate,
                        parentTrace.assumeblocked);
}


/* *************************** */
/* EXTEND TRACE                */
/* *************************** */

ZExplorer::TraceExtension
ZExplorer::extendTrace(std::vector<ZEvent>&& tr)
{
  clock_t init = std::clock();
  ZBuilderTSO TB(originalTB.config, originalTB.M, std::move(tr));
  auto traceExtension = TraceExtension(TB.extendGivenTrace());
  time_interpreter += (double)(clock() - init)/CLOCKS_PER_SEC;
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
/* RESPECTS ANNOTATION         */
/* *************************** */

bool ZExplorer::respectsAnnotation
(const std::vector<ZEvent>& trace,
 const ZAnnotation& annotation,
 const ZTrace& parentTrace) const
{
  for (const ZEvent& evref: trace) {
    const ZEvent *ev = &evref;
    if (isRead(ev) && annotation.defines(ev->threadID(), ev->eventID())) {
      const ZObs& obs = annotation.getObs(ev);
      const ZEvent *obsB = (obs.thr == INT_MAX)
        ? parentTrace.graph.getBasis().initial() : parentTrace.graph.getBasis().getEvent(obs);
      const ZEvent *obsM = (isInitial(obsB))
        ? parentTrace.graph.getBasis().initial() : &(trace.at(obsB->write_other_trace_id));
      assert(obsB->value == obsM->value);
      assert(isInitial(obsB) || (isWriteB(obsB) && isWriteM(obsM) &&
                                 sameMl(obsB, obsM) && sameMl(ev, obsB)));
      const ZEvent *realObservation = (ev->observed_trace_id == -1)
        ? parentTrace.graph.getBasis().initial() : &(trace.at(ev->observed_trace_id));

      if (realObservation != obsB && realObservation != obsM) {
        parentTrace.dump();
        dumpTrace(trace);
        annotation.dump();
        llvm::errs() << "This read         :::  ";
        ev->dump();
        llvm::errs() << "Should observeB   :::  ";
        obsB->dump();
        llvm::errs() << "Should observeM   :::  ";
        obsM->dump();
        llvm::errs() << "Actually observed :::  ";
        realObservation->dump();
        return false;
      }
    }
  }

  return true;
}

/*

  std::vector<unsigned> processLengths = current->graph.getProcessLengths();
  auto negativeWriteMazBranch = ZAnnotationNeg(current->negative);
  negativeWriteMazBranch.update(nd, processLengths);

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
*/
