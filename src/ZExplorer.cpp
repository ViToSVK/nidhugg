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
  : originalTB(&tb), tso(true)
{
  if (tb.someThreadAssumeBlocked)
    assume_blocked_thread = 1;

  if (tb.somethingToAnnotate.empty())
    executed_traces_full = 1;
  else
    initial = new ZTrace(std::move(tb.prefix), tb.star_root_index,
                         tb.someThreadAssumeBlocked, tso);
}


ZExplorer::ZExplorer(ZBuilderPSO& tb)
  : originalTB(&tb), tso(false)
{
  if (tb.someThreadAssumeBlocked)
    assume_blocked_thread = 1;

  if (tb.somethingToAnnotate.empty())
    executed_traces_full = 1;
  else
    initial = new ZTrace(std::move(tb.prefix), tb.star_root_index,
                         tb.someThreadAssumeBlocked, tso);
}


void ZExplorer::print_stats() const
{
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
  std::cout << "Closure succ -- no added edge:     " << closure_no_edge << "\n";
  std::cout << "Closure succ -- total added edges: " << closure_edges << "\n";
  std::cout << "Closure succ -- total iterations:  " << closure_iter << "\n";
  std::cout << std::setprecision(2) << std::fixed;
  std::cout << "Time spent on copying:             " << time_copy << "\n";
  std::cout << "Time spent on linearization:       " << time_linearization << "\n";
  std::cout << "Time spent on interpreting:        " << time_interpreter << "\n";
  std::cout << "Time spent on chronological:       " << time_chrono << "\n";
  std::cout << "Time spent on closure:             " << time_closure << "\n";
  std::cout << "Time spent on closure-succ-noedge: " << time_closure_no_edge << "\n";
  std::cout << "\n" << std::scientific;

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
  if (info) {
    //dumpTrace(annTrace.trace);
    annTrace.graph.dump();
    annTrace.annotation.dump();
    //llvm::errs() << "-------------------------\n\n";
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
        assert(originalTB && originalTB->error_trace);
        return error;
      }
    }
    else {
      assert(isLock(ev));
      bool error = mutateLock(annTrace, ev);
      if (error) {
        assert(originalTB && originalTB->error_trace);
        return error;
      }
    }
    annTrace.negative.update
      (ev, annTrace.graph.getBasis().real_sizes_minus_one());
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
    if (info) {
      llvm::errs() << "Read to mutate:\n";
      read->dump();
      llvm::errs() << "All candidates forbidden\n";
      annTrace.annotation.dump();
      annTrace.negative.dump();
    }
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
      annTrace.negative.dump();
    }

    auto init = std::clock();
    ZClosure preClosure(mutatedAnnotation, mutatedPO);
    preClosure.preClose(read, observation);
    time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

    bool mutationFollowsCurrentTrace = false;
    if (observation.thr == INT_MAX) {
      mutationFollowsCurrentTrace = (read->value == 0);
      // DC (read->observed_trace_id == -1);
    } else {
      const ZEvent *obsB = annTrace.getEvent(observation);
      assert(isWriteB(obsB) && isWriteM(obsB->write_other_ptr) &&
             sameMl(read, obsB) && sameMl(read, obsB->write_other_ptr));
      mutationFollowsCurrentTrace = (read->value == obsB->value);
      // DC (read->observed_trace_id == obsB->traceID() ||
      // DC  read->observed_trace_id == obsB->write_other_ptr->traceID());
    }

    bool error = chronological
      (annTrace, read, std::move(mutatedAnnotation),
       std::move(mutatedPO), mutationFollowsCurrentTrace);

    if (error) {
      assert(originalTB && originalTB->error_trace);
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
      annTrace.negative.dump();
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
    annTrace.negative.dump();
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
       std::move(mutatedPO), mutationFollowsCurrentTrace,
       nullptr);
  }

  auto toOrder = annTrace.graph.chronoOrderPairs
    ((isRead(readLock) && !annTrace.isRoot(readLock) &&
      annTrace.isRoot(mutatedAnnotation.getObs(readLock).thr))
     ? readLock : nullptr, mutatedAnnotation);

  if (toOrder.empty()) {
    // No pairs to order
    return chronoReads
      (annTrace, readLock, std::move(mutatedAnnotation),
       std::move(mutatedPO), mutationFollowsCurrentTrace,
       nullptr);
  }

  // First - partial order
  // Second - check toOrder pairs starting from this one
  std::list<std::pair<ZPartialOrder, unsigned>> worklist;
  worklist.emplace_front(std::move(mutatedPO), 0);
  assert(mutatedPO.empty() && !worklist.empty() &&
         !worklist.front().first.empty());

  bool sawFullTrace = false;
  while (!worklist.empty()) {
    if (sawFullTrace) {
      // This mutation leads to a full trace without an assertion violation
      // all the other successful chrono orderings generated here
      // would produce traces with the same events, hence no need to check them
      break;
    }

    // Recursively process one chronoPO branch
    auto init = std::clock();
    std::pair<ZPartialOrder, unsigned> current
      (std::move(worklist.front()));
    assert(!current.first.empty() && !worklist.empty() &&
           worklist.front().first.empty());
    worklist.pop_front();

    while (current.second < toOrder.size()) {
      auto ev1 = toOrder[current.second].first;
      auto ev2 = toOrder[current.second].second;
      // Order pair ev1-ev2
      if (!current.first.areOrdered(ev1, ev2)) {
        // *one-read*-or-*both-mws-observable-in-po* condition
        if (isRead(ev1) || isRead(ev2) ||
            (annTrace.graph.isObservable(ev1, current.first) &&
             annTrace.graph.isObservable(ev2, current.first))) {
          // Create two cases with these orderings
          worklist.emplace_front(ZPartialOrder(current.first), // copy
                                 current.second + 1);
          // Handle current: ev1 -> ev2
          current.first.addEdge(ev1, ev2);
          // Handle otherorder: ev2 -> ev1
          worklist.front().first.addEdge(ev2, ev1);
        }
      }
      current.second++;
    }
    assert(current.second == toOrder.size());
    #ifndef NDEBUG
    for (auto& pair : toOrder)
      assert(current.first.areOrdered(pair.first, pair.second) ||
             (!isRead(pair.first) && !isRead(pair.second) &&
              (!annTrace.graph.isObservable(pair.first, current.first) ||
               !annTrace.graph.isObservable(pair.second, current.first))));
      // *ordered*-or-(*no-read*-and-*some-mw-unobservable-in-po*)
    #endif
    time_chrono += (double)(clock() - init)/CLOCKS_PER_SEC;

    bool error = chronoReads
      (annTrace, readLock,
       worklist.empty()
       ? std::move(mutatedAnnotation) // move
       : ZAnnotation(mutatedAnnotation), // copy
       std::move(current.first), mutationFollowsCurrentTrace,
       &sawFullTrace);

    if (error) {
      assert(originalTB && originalTB->error_trace);
      return error;
    }
  }
  worklist.clear();

  return false;
}


bool ZExplorer::chronoReads
(const ZTrace& annTrace, const ZEvent *readLock,
 ZAnnotation&& mutatedAnnotation, ZPartialOrder&& mutatedPO,
 bool mutationFollowsCurrentTrace, bool *fullTraceAfterChrono)
{
  // This PO has all conflicting leaf memory-write pairs ordered
  // Here we have to ensure that all leaf reads with leaf observation
  // happen before the leaf memory-write right after the observation-memory-write
  // (no matter if the observation is local-leaf or remote-leaf)

  auto init = std::clock();

  std::list<const ZEvent *> leafReadsLeafObs;
  if (isRead(readLock) && !annTrace.isRoot(readLock) &&
      !annTrace.isRoot(mutatedAnnotation.getObs(readLock).thr))
    leafReadsLeafObs.push_back(readLock);
  for (const auto& thr_list : annTrace.graph.getCache().chronoAnnR)
    for (const auto& leaf : thr_list.second) {
      assert(isRead(leaf) && mutatedAnnotation.defines(leaf) &&
             !annTrace.isRoot(leaf) && leaf->threadID() == thr_list.first);
      if (!annTrace.isRoot(mutatedAnnotation.getObs(leaf).thr))
        leafReadsLeafObs.push_back(leaf);
    }

  for (const auto& leaf : leafReadsLeafObs) {
    const ZEvent *anBuf = mutatedAnnotation.getObs(leaf).thr != INT_MAX
      ? annTrace.getEvent(mutatedAnnotation.getObs(leaf)) : nullptr;
    assert(!anBuf || isWriteB(anBuf));
    assert(!anBuf || !annTrace.isRoot(anBuf)); // leaf observation
    const ZEvent *anMem = anBuf ? anBuf->write_other_ptr : nullptr;
    assert(!anMem ||
           (isWriteM(anMem) && sameMl(anMem, leaf) &&
            (leaf->threadID() == anMem->threadID() ||
             mutatedPO.hasEdge(anMem, leaf)))); // PreClosure

    // Find the earliest conflicting read-remote leaf memory-write
    // happening after the observation-memory-write
    const ZEvent *nextMem = nullptr;
    for (const auto& thr_list : annTrace.graph.getCache().chrono) { // leaf
      if (thr_list.first != leaf->threadID()) { // read-remote
        for (const auto& mwEv : thr_list.second) {

          assert(isWriteM(mwEv) && mwEv->threadID() == thr_list.first);
          if (mwEv != anMem && sameMl(mwEv, leaf)) { // conflicting
            assert(!anMem || mutatedPO.areOrdered(anMem, mwEv) ||
                   !annTrace.graph.isObservable(anMem, mutatedPO) ||
                   !annTrace.graph.isObservable(mwEv, mutatedPO));
            if (!anMem || mutatedPO.hasEdge(anMem, mwEv)) {
              // all conflicting others in this list happen after mwEv
              // so only mwEv can be a candidate; break
              if (annTrace.graph.isObservable(mwEv, mutatedPO)) {
                if (!nextMem)
                  nextMem = mwEv;
                else {
                  assert(mutatedPO.areOrdered(nextMem, mwEv));
                  if (mutatedPO.hasEdge(mwEv, nextMem))
                    nextMem = mwEv;
                }
              }
              break;
            }
          }

        } // break from this loop
      }
    }
    // the read has to happen before nextMem
    if (nextMem && mutatedPO.hasEdge(nextMem, leaf)) {
      time_chrono += (double)(clock() - init)/CLOCKS_PER_SEC;
      return false;
    }
    else if (nextMem && !mutatedPO.hasEdge(leaf, nextMem))
      mutatedPO.addEdge(leaf, nextMem);
  }

  time_chrono += (double)(clock() - init)/CLOCKS_PER_SEC;
  ++leaf_chrono_pos;
  return closePO
    (annTrace, readLock, std::move(mutatedAnnotation),
     std::move(mutatedPO), mutationFollowsCurrentTrace,
     fullTraceAfterChrono);
}


/* *************************** */
/* CLOSE PO                    */
/* *************************** */

bool ZExplorer::closePO
(const ZTrace& annTrace, const ZEvent *readLock,
 ZAnnotation&& mutatedAnnotation, ZPartialOrder&& mutatedPO,
 bool mutationFollowsCurrentTrace, bool *fullTraceAfterChrono)
{
  auto init = std::clock();
  ZClosure closure(mutatedAnnotation, mutatedPO);
  bool closed = closure.close
    (isLock(readLock) ? nullptr : readLock);
  double time = (double)(clock() - init)/CLOCKS_PER_SEC;
  time_closure += time;

  if (!closed) {
    //if (info) llvm::errs() << "Closure failed\n\n-----\n\n";
    ++closure_failed;
    return false;
  }

  //if (info) llvm::errs() << "Closure succeeded\n\n";
  ++closure_succeeded;
  if (closure.added_edges == 0) {
    closure_no_edge++;
    time_closure_no_edge += time;
  }
  closure_edges += closure.added_edges;
  closure_iter += closure.iterations;

  return extendAndRecur
    (annTrace, std::move(mutatedAnnotation),
     std::move(mutatedPO), mutationFollowsCurrentTrace,
     fullTraceAfterChrono);
}


/* *************************** */
/* EXTEND AND RECUR            */
/* *************************** */

bool ZExplorer::extendAndRecur
(const ZTrace& parentTrace, ZAnnotation&& mutatedAnnotation,
 ZPartialOrder&& mutatedPO, bool mutationFollowsCurrentTrace,
 bool *fullTraceAfterChrono)
{
  TraceExtension mutatedTrace;
  if (mutationFollowsCurrentTrace)
    mutatedTrace = reuseTrace(parentTrace, mutatedAnnotation);
  else {
    clock_t init = std::clock();
    ZLinearization linearizer(mutatedAnnotation, mutatedPO);
    auto linear = tso ? linearizer.linearizeTSO() : linearizer.linearizePSO();
    time_linearization += (double)(clock() - init)/CLOCKS_PER_SEC;
    assert(linearizationRespectsAnn(linear, mutatedAnnotation, mutatedPO, parentTrace));

    mutatedTrace = extendTrace(std::move(linear));
    // if (info) dumpTrace(mutatedTrace.trace);
    assert(respectsAnnotation(mutatedTrace.trace, mutatedAnnotation, mutatedPO, parentTrace));
  }

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
    if (info) {
      llvm::errs() << "FULL\n";
      //mutatedAnnotation.dump();
    }
    // Note in explorer that mutation produces max-trace
    if (fullTraceAfterChrono)
      *fullTraceAfterChrono = true;
    return false;
  }

  clock_t init = std::clock();
  assert(!mutatedTrace.empty() && !mutatedPO.empty());
  ZTrace mutatedZTrace
    (parentTrace,
     std::move(mutatedTrace.trace),
     std::move(mutatedAnnotation),
     std::move(mutatedPO),
     mutatedTrace.hasAssumeBlockedThread);
  assert(mutatedTrace.empty() && mutatedPO.empty() &&
         mutatedAnnotation.empty());
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
  assert(originalTB);
  if (tso) {
    clock_t init = std::clock();
    ZBuilderTSO TB(*(originalTB->config), originalTB->M, std::move(tr));
    auto traceExtension = TraceExtension(TB.extendGivenTrace());
    time_interpreter += (double)(clock() - init)/CLOCKS_PER_SEC;
    interpreter_used++;

    if (TB.has_error()) {
      // ERROR FOUND
      originalTB->error_trace = TB.get_trace();
      traceExtension.hasError = true;
    }

    if (TB.someThreadAssumeBlocked)
      traceExtension.hasAssumeBlockedThread = true;

    return traceExtension;
  }

  assert(!tso);
  clock_t init = std::clock();
  ZBuilderPSO TB(*(originalTB->config), originalTB->M, std::move(tr));
  auto traceExtension = TraceExtension(TB.extendGivenTrace());
  time_interpreter += (double)(clock() - init)/CLOCKS_PER_SEC;
  interpreter_used++;

  if (TB.has_error()) {
    // ERROR FOUND
    originalTB->error_trace = TB.get_trace();
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
 const ZPartialOrder& mutatedPO,
 const ZTrace& parentTrace) const
{
  assert(!trace.empty());
  std::unordered_map<ZObs, unsigned> bw_pos;
  for (const ZEvent& evref: trace) {
    const ZEvent *ev = &evref;
    // ev->threadID() not set, have to get it
    unsigned thrid = parentTrace.graph.getBasis().getThreadIDnoAdd(ev);
    ev->_thread_id = thrid;
    if (isWriteB(ev) && thrid != 1337)
      bw_pos.emplace(ZObs(thrid, ev->eventID()), ev->traceID());
  }
  for (const ZEvent& evref : trace) {
    const ZEvent *ev = &evref;
    if (isRead(ev)) {
      if (annotation.defines(ev->threadID(), ev->eventID())) {
        const ZObs& obs = annotation.getObs(ev->threadID(), ev->eventID());
        const ZEvent *obsB = parentTrace.graph.getBasis().initial();
        if (obs.thr != INT_MAX) {
          assert(bw_pos.count(obs));
          obsB = &(trace.at(bw_pos.at(obs)));
        }
        const ZEvent *obsM = (isInitial(obsB))
          ? parentTrace.graph.getBasis().initial() : &(trace.at(obsB->write_other_trace_id));
        assert(obsB->value == obsM->value);
        assert(isInitial(obsB) || (isWriteB(obsB) && isWriteM(obsM) &&
                                   sameMl(obsB, obsM) && sameMl(ev, obsB)));
        const ZEvent *realObservation = (ev->observed_trace_id == -1)
          ? parentTrace.graph.getBasis().initial() : &(trace.at(ev->observed_trace_id));
        if (realObservation != obsB && realObservation != obsM) {
          parentTrace.dump();
          llvm::errs() << "Closed annotated partial order\n";
          mutatedPO.dump();
          llvm::errs() << "Extension\n";
          dumpTrace(trace);
          llvm::errs() << "Full annotation that should be respected\n";
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
  }
  // Unset thread-id, let Graph take care of it
  for (const ZEvent& evref : trace)
    evref._thread_id = 1337;
  return true;
}


bool ZExplorer::linearizationRespectsAnn
(const std::vector<ZEvent>& trace,
 const ZAnnotation& annotation,
 const ZPartialOrder& mutatedPO,
 const ZTrace& parentTrace) const
{
  assert(!trace.empty());
  // Trace index of the last write happening in the given location
  std::unordered_map<SymAddrSize, int> lastWrite;
  // ThreadID -> ML -> Writes in store queue of thr for ml
  std::unordered_map
    <int, std::unordered_map
     <SymAddrSize, std::list<int>>> storeQueue;
  // Real observation in trace
  std::unordered_map<int, int> realObs;
  // BufferWrite -> pos
  std::unordered_map<ZObs, unsigned> bw_pos;
  // BufferWrite -> MemoryWrite
  std::unordered_map<int, int> buf_mem;

  assert(!trace.empty());
  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = &(trace.at(i));
    ev->_thread_id = parentTrace.graph.getBasis().getThreadIDnoAdd(ev);
    assert(ev->threadID() < 1337);

    if (isWriteB(ev)) {
      bw_pos.emplace(ZObs(ev->threadID(), ev->eventID()), i);
      if (!storeQueue.count(ev->threadID()))
        storeQueue.emplace
          (ev->threadID(), std::unordered_map<SymAddrSize, std::list<int>>());
      if (!storeQueue.at(ev->threadID()).count(ev->ml))
        storeQueue.at(ev->threadID()).emplace
          (ev->ml, std::list<int>());
      storeQueue.at(ev->threadID()).at(ev->ml).push_back(i);
    }

    if (isWriteM(ev)) {
      lastWrite[ev->ml] = i;
      assert(!storeQueue.at(ev->threadID()).at(ev->ml).empty());
      buf_mem.emplace(storeQueue.at(ev->threadID()).at(ev->ml).front(), i);
      storeQueue.at(ev->threadID()).at(ev->ml).pop_front();
    }

    if (isRead(ev)) {
      if (storeQueue.count(ev->threadID()) &&
          storeQueue.at(ev->threadID()).count(ev->ml) &&
          !storeQueue.at(ev->threadID()).at(ev->ml).empty())
        realObs.emplace(i, storeQueue.at(ev->threadID()).at(ev->ml).back());
      else if (lastWrite.count(ev->ml))
        realObs.emplace(i, lastWrite.at(ev->ml));
      else
        realObs.emplace(i, -1);
    }
  }

  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = &(trace.at(i));
    if (isRead(ev)) {
      assert(annotation.defines(ev->threadID(), ev->eventID()));
      const ZObs& obs = annotation.getObs(ev->threadID(), ev->eventID());
      const ZEvent *obsB = parentTrace.graph.getBasis().initial();
      if (obs.thr != INT_MAX) {
        assert(bw_pos.count(obs));
        obsB = &(trace.at(bw_pos.at(obs)));
        assert(isWriteB(obsB));
      }
      const ZEvent *obsM = (isInitial(obsB))
        ? parentTrace.graph.getBasis().initial() : &(trace.at(buf_mem.at(bw_pos.at(obs))));;
      assert(obsB->value == obsM->value);
      assert(isInitial(obsB) || (isWriteB(obsB) && isWriteM(obsM) &&
                                 sameMl(obsB, obsM) && sameMl(ev, obsB)));
      const ZEvent *realObservation = (realObs.at(i) == -1)
        ? parentTrace.graph.getBasis().initial() : &(trace.at(realObs.at(i)));
      if (realObservation != obsB && realObservation != obsM) {
        parentTrace.dump();
        llvm::errs() << "Closed annotated partial order\n";
        mutatedPO.dump();
        llvm::errs() << "Linearization\n";
        dumpTrace(trace);
        llvm::errs() << "Full annotation that should be respected\n";
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

  // Unset thread-id, let Graph take care of it
  for (const ZEvent& evref : trace)
    evref._thread_id = 1337;

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
