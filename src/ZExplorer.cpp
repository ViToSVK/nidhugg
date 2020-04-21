/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2020 Viktor Toman
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

static const bool DEBUG = false;
#include "ZDebug.h"


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


ZExplorer::ZExplorer(ZBuilderSC& tb)
  : originalTB(&tb), tso(true), sc_flag(true)
{
  if (tb.someThreadAssumeBlocked)
    assume_blocked_thread = 1;

  if (tb.somethingToAnnotate.empty())
    executed_traces_full = 1;
  else
    initial = new ZTrace(std::move(tb.prefix), tb.someThreadAssumeBlocked, tso);
}


ZExplorer::ZExplorer(ZBuilderTSO& tb)
  : originalTB(&tb), tso(true), sc_flag(false)
{
  if (tb.someThreadAssumeBlocked)
    assume_blocked_thread = 1;

  if (tb.somethingToAnnotate.empty())
    executed_traces_full = 1;
  else
    initial = new ZTrace(std::move(tb.prefix), tb.someThreadAssumeBlocked, tso);
}


ZExplorer::ZExplorer(ZBuilderPSO& tb)
  : originalTB(&tb), tso(false), sc_flag(false)
{
  if (tb.someThreadAssumeBlocked)
    assume_blocked_thread = 1;

  if (tb.somethingToAnnotate.empty())
    executed_traces_full = 1;
  else
    initial = new ZTrace(std::move(tb.prefix), tb.someThreadAssumeBlocked, tso);
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
  std::cout << "Closure failed:                    " << closure_failed << "\n";
  std::cout << "Closure succeeded:                 " << closure_succeeded << "\n";
  std::cout << "Closure succ -- no added edge:     " << closure_no_edge << "\n";
  std::cout << "Closure succ -- total added edges: " << closure_edges << "\n";
  std::cout << "Closure succ -- total iterations:  " << closure_iter << "\n";
  std::cout << std::setprecision(2) << std::fixed;
  std::cout << "Time spent on copying:             " << time_copy << "\n";
  std::cout << "Time spent on linearization:       " << time_linearization << "\n";
  std::cout << "Time spent on interpreting:        " << time_interpreter << "\n";
  std::cout << "Time spent on closure:             " << time_closure << "\n";
  std::cout << "Time spent on closure-succ-noedge: " << time_closure_no_edge << "\n";
  std::cout << "Linearization branching factor:    " << (double)total_children/total_parents << "\n";
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
  #ifndef NDEBUG
    std::cout << "RUNNING DEBUG VERSION\n";
  #endif

  if (initial)
    return exploreRec(*initial);
  return false;
}


bool ZExplorer::exploreRec(ZTrace& annTrace)
{
  start_err("exploreRec...");
  if (info) {
    //dumpTrace(annTrace.trace);
    annTrace.graph.dump();
    annTrace.annotation.dump();
    annTrace.negative.dump();
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
        end_err("?a");
        return error;
      }
    }
    else {
      assert(isLock(ev));
      bool error = mutateLock(annTrace, ev);
      if (error) {
        assert(originalTB && originalTB->error_trace);
        end_err("?b");
        return error;
      }
    }
    annTrace.negative.update
      (ev, annTrace.graph.real_sizes_minus_one());
  }

  // Done with this recursion node
  if (annTrace.deadlocked) {
    // Not full trace ending in a deadlock
    // We count it as full
    ++executed_traces_full;
    ++executed_traces_full_deadlock;
  }

  end_err("0");
  return false;
}


/* *************************** */
/* MUTATE READ                 */
/* *************************** */

bool ZExplorer::mutateRead(const ZTrace& annTrace, const ZEvent *read)
{
  start_err("mutateRead...");
  assert(isRead(read));
  if (info) {
      llvm::errs() << "Read to mutate:\n";
      read->dump();
  }

  auto obsCandidates = annTrace.getObsCandidates(read);
  if (obsCandidates.empty()) {
    // All candidates are ruled out by negative annotation
    ++no_mut_choices;
    if (info) {
      llvm::errs() << "All candidates forbidden\n";
    }
    end_err("0a");
    return false;
  }

  for (auto& observation : obsCandidates) {
    ++mutations_considered;

    ZAnnotation mutatedAnnotation(annTrace.annotation);
    mutatedAnnotation.add(read, observation);
    assert(mutatedAnnotation.size() == annTrace.annotation.size() + 1);
    auto mutatedPO = annTrace.graph.copyPO();

    if (info) {
      llvm::errs() << "Observation:\n" << observation.to_string() << "\n";
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
      // DC (read->observed_trace_id == obsB->trace_id() ||
      // DC  read->observed_trace_id == obsB->write_other_ptr->trace_id());
    }

    bool error = closePO
      (annTrace, read, std::move(mutatedAnnotation),
       std::move(mutatedPO), mutationFollowsCurrentTrace);

    if (error) {
      assert(originalTB && originalTB->error_trace);
      end_err("1");
      return error;
    }
  }

  end_err("0b");
  return false;
}


/* *************************** */
/* MUTATE LOCK                 */
/* *************************** */

bool ZExplorer::mutateLock(const ZTrace& annTrace, const ZEvent *lock)
{
  start_err("mutateLock...");
  assert(isLock(lock));

  if (info) {
    llvm::errs() << "Lock to mutate:\n";
    lock->dump();
  }

  if (!annTrace.annotation.locationHasSomeLock(lock)) {
    // This lock hasn't been touched before
    annTrace.deadlocked = false;
    if (annTrace.negative.forbidsInitialEvent(lock)) {
      // Negative annotation forbids initial unlock
      end_err("0a");
      if (info) {
        llvm::errs() << "Negative annotation forbids initial unlock\n";
      }
      return false;
    }

    // Trivially realizable
    ++mutations_considered;
    ZAnnotation mutatedAnnotation(annTrace.annotation);
    mutatedAnnotation.setLastLock(lock);
    auto mutatedPO = annTrace.graph.copyPO();

    if (info) {
      llvm::errs() << "Not locked before\n";
    }

    bool mutationFollowsCurrentTrace =
      (lock->observed_trace_id == -1);

    auto res = closePO
      (annTrace, lock, std::move(mutatedAnnotation),
       std::move(mutatedPO), mutationFollowsCurrentTrace);
    end_err("?a");
    return res;
  }

  // The lock has been touched before
  auto lastLockObs = annTrace.annotation.getLastLock(lock);
  auto lastUnlock = annTrace.graph.getUnlockOfThisLock(lastLockObs);

  if (!lastUnlock) {
    // This lock is currently locked
    // Trivially unrealizable
    end_err("0b");
    if (info) {
      llvm::errs() << "Currently locked\n";
    }
    return false;
  }

  // This lock is currently unlocked by lastUnlock
  assert(lastUnlock && isUnlock(lastUnlock) &&
         sameMl(lock, lastUnlock));
  annTrace.deadlocked = false;

  if (annTrace.negative.forbids(lock, lastUnlock)) {
    // Negative annotation forbids this unlock
    end_err("0c");
    return false;
  }

  // Realizable
  ++mutations_considered;
  ZAnnotation mutatedAnnotation(annTrace.annotation);
  mutatedAnnotation.setLastLock(lock);
  auto mutatedPO = annTrace.graph.copyPO();

  if (info) {
    llvm::errs() << "Currently unlocked by:\n";
    lastUnlock->dump();
  }

  auto init = std::clock();
  ZClosure preClosure(mutatedAnnotation, mutatedPO);
  preClosure.preClose(lock, lastUnlock);
  time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

  bool mutationFollowsCurrentTrace =
    (lock->observed_trace_id == lastUnlock->trace_id());

  end_err("?b");
  return closePO
    (annTrace, lock, std::move(mutatedAnnotation),
     std::move(mutatedPO), mutationFollowsCurrentTrace);
}


/* *************************** */
/* CLOSE PO                    */
/* *************************** */

bool ZExplorer::closePO
(const ZTrace& annTrace, const ZEvent *readLock,
 ZAnnotation&& mutatedAnnotation, ZPartialOrder&& mutatedPO,
 bool mutationFollowsCurrentTrace)
{
  start_err("closePO...");
  auto init = std::clock();
  ZClosure closure(mutatedAnnotation, mutatedPO);
  bool closed = closure.close
    (isLock(readLock) ? nullptr : readLock);
  double time = (double)(clock() - init)/CLOCKS_PER_SEC;
  time_closure += time;

  if (!closed) {
    //if (info) llvm::errs() << "Closure failed\n\n-----\n\n";
    ++closure_failed;
    end_err("0");
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

  auto res = extendAndRecur
    (annTrace, std::move(mutatedAnnotation),
     std::move(mutatedPO), mutationFollowsCurrentTrace);
  end_err("?");
  return res;
}


/* *************************** */
/* EXTEND AND RECUR            */
/* *************************** */

bool ZExplorer::extendAndRecur
(const ZTrace& parentTrace, ZAnnotation&& mutatedAnnotation,
 ZPartialOrder&& mutatedPO, bool mutationFollowsCurrentTrace)
{
  start_err("extendAndRecur...");
  TraceExtension mutatedTrace;
  if (mutationFollowsCurrentTrace) {
    mutatedTrace = reuseTrace(parentTrace, mutatedAnnotation);
  }
  else {
    clock_t init = std::clock();
    ZLinearization linearizer(mutatedAnnotation, mutatedPO, parentTrace.trace);
    auto linear = tso ? linearizer.linearizeTSO() : linearizer.linearizePSO();
    time_linearization += (double)(clock() - init)/CLOCKS_PER_SEC;
    total_parents += linearizer.num_parents;
    total_children += linearizer.num_children;
    if (linear.empty()) {
      end_err("0a");
      return false;
    }
    assert(linearizationRespectsAnn(linear, mutatedAnnotation, mutatedPO, parentTrace));
    // if (info) dumpTrace(linear);
    mutatedTrace = extendTrace(std::move(linear));
    // if (info) dumpTrace(mutatedTrace.trace);
    assert(respectsAnnotation(mutatedTrace.trace, mutatedAnnotation, mutatedPO, parentTrace));
  }

  executed_traces++;
  if (mutatedTrace.hasError) {
    assert(!mutationFollowsCurrentTrace);
    end_err("1");
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
    end_err("0b");
    return false;
  }

  clock_t init = std::clock();
  assert(!mutatedTrace.empty() && !mutatedPO.empty());
  err_msg("creating extension");
  ZTrace mutatedZTrace
    (parentTrace,
     std::move(mutatedTrace.trace),
     std::move(mutatedAnnotation),
     std::move(mutatedPO),
     mutatedTrace.hasAssumeBlockedThread);
  err_msg("created extension");
  assert(mutatedTrace.empty() && mutatedPO.empty() &&
         mutatedAnnotation.empty());
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;

  auto res = exploreRec(mutatedZTrace);
  end_err("?");
  return res;
}


/* *************************** */
/* REUSE TRACE                 */
/* *************************** */

ZExplorer::TraceExtension
ZExplorer::reuseTrace
(const ZTrace& parentTrace, const ZAnnotation& mutatedAnnotation)
{
  start_err("reuseTrace..");
  std::vector<ZEvent> tr;
  tr.reserve(parentTrace.trace.size());
  bool somethingToAnnotate = false;

  for (const ZEvent& ev : parentTrace.trace) {
    if (!somethingToAnnotate) {
      if (isRead(ev) && (!mutatedAnnotation.defines(&ev)))
        somethingToAnnotate = true;
      else if (isLock(ev)) {
        if (!parentTrace.graph.hasEvent(&ev))
          somethingToAnnotate = true;
        else if (ev.event_id() == // last in its thraux
                 (parentTrace.graph)(ev.thread_id(), ev.aux_id()).size() - 1
                 && !mutatedAnnotation.isLastLock(&ev))
          somethingToAnnotate = true;
      }
    }

    tr.push_back(ev.copy(tr.size(), false));
  }

  auto res = TraceExtension(std::move(tr), somethingToAnnotate,
                            parentTrace.assumeblocked);
  end_err("?");
  return res;
}


/* *************************** */
/* EXTEND TRACE                */
/* *************************** */

ZExplorer::TraceExtension
ZExplorer::extendTrace(std::vector<ZEvent>&& tr)
{
  start_err("extendTrace...");
  assert(originalTB);
  if (tso && sc_flag) {
    clock_t init = std::clock();
    ZBuilderSC TB(*(originalTB->config), originalTB->M, std::move(tr));
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

    end_err("?a");
    return traceExtension;
  }
  else if (tso && !sc_flag) {
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

    end_err("?a");
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

  end_err("?b");
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
    // ev->thread_id() not set, have to get it
    unsigned thrid = parentTrace.graph.getThreadIDnoAdd(ev);
    ev->_thread_id = thrid;
    if (isWriteB(ev) && thrid != 1337)
      bw_pos.emplace(ZObs(thrid, ev->event_id()), ev->trace_id());
  }
  for (const ZEvent& evref : trace) {
    const ZEvent *ev = &evref;
    if (isRead(ev)) {
      if (annotation.defines(ev->thread_id(), ev->event_id())) {
        const ZObs& obs = annotation.getObs(ev->thread_id(), ev->event_id());
        const ZEvent *obsB = parentTrace.graph.initial();
        if (obs.thr != INT_MAX) {
          assert(bw_pos.count(obs));
          obsB = &(trace.at(bw_pos.at(obs)));
        }
        const ZEvent *obsM = (isInitial(obsB))
          ? parentTrace.graph.initial() : &(trace.at(obsB->write_other_trace_id));
        assert(obsB->value == obsM->value);
        assert(isInitial(obsB) || (isWriteB(obsB) && isWriteM(obsM) &&
                                   sameMl(obsB, obsM) && sameMl(ev, obsB)));
        const ZEvent *realObservation = (ev->observed_trace_id == -1)
          ? parentTrace.graph.initial() : &(trace.at(ev->observed_trace_id));
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
  // Trace index of the last lock + whether given ML is currently locked
  std::unordered_map<SymAddrSize, int> lastLock;
  std::unordered_set<SymAddrSize> locked;
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

  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = &(trace.at(i));
    ev->_thread_id = parentTrace.graph.getThreadIDnoAdd(ev);
    assert(ev->thread_id() < 1337);

    if (isWriteB(ev)) {
      bw_pos.emplace(ZObs(ev->thread_id(), ev->event_id()), i);
      if (!storeQueue.count(ev->thread_id()))
        storeQueue.emplace
          (ev->thread_id(), std::unordered_map<SymAddrSize, std::list<int>>());
      if (!storeQueue.at(ev->thread_id()).count(ev->ml))
        storeQueue.at(ev->thread_id()).emplace
          (ev->ml, std::list<int>());
      storeQueue.at(ev->thread_id()).at(ev->ml).push_back(i);
    }

    if (isWriteM(ev)) {
      lastWrite[ev->ml] = i;
      assert(storeQueue.count(ev->thread_id()) && "Maybe writeM before writeB?");
      assert(storeQueue.at(ev->thread_id()).count(ev->ml) && "Maybe writeM before writeB?");
      assert(!storeQueue.at(ev->thread_id()).at(ev->ml).empty());
      buf_mem.emplace(storeQueue.at(ev->thread_id()).at(ev->ml).front(), i);
      storeQueue.at(ev->thread_id()).at(ev->ml).pop_front();
    }

    if (isRead(ev)) {
      if (storeQueue.count(ev->thread_id()) &&
          storeQueue.at(ev->thread_id()).count(ev->ml) &&
          !storeQueue.at(ev->thread_id()).at(ev->ml).empty())
        realObs.emplace(i, storeQueue.at(ev->thread_id()).at(ev->ml).back());
      else if (lastWrite.count(ev->ml))
        realObs.emplace(i, lastWrite.at(ev->ml));
      else
        realObs.emplace(i, -1);
    }

    if (isLock(ev)) {
      assert(!locked.count(ev->ml));
      locked.insert(ev->ml);
      lastLock[ev->ml] = i;
    }

    if (isUnlock(ev)) {
      assert(locked.count(ev->ml));
      locked.erase(ev->ml);
    }
  }

  // Check whether last locks are consistent with annotation
  for (auto entry : lastLock) {
    const ZEvent *ev = &(trace.at(entry.second));
    assert(annotation.isLastLock(ev));
  }

  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = &(trace.at(i));
    if (isRead(ev)) {
      assert(annotation.defines(ev->thread_id(), ev->event_id()));
      const ZObs& obs = annotation.getObs(ev->thread_id(), ev->event_id());
      const ZEvent *obsB = parentTrace.graph.initial();
      if (obs.thr != INT_MAX) {
        assert(bw_pos.count(obs));
        obsB = &(trace.at(bw_pos.at(obs)));
        assert(isWriteB(obsB));
      }
      const ZEvent *obsM = (isInitial(obsB))
        ? parentTrace.graph.initial() : &(trace.at(buf_mem.at(bw_pos.at(obs))));;
      assert(obsB->value == obsM->value);
      assert(isInitial(obsB) || (isWriteB(obsB) && isWriteM(obsM) &&
                                 sameMl(obsB, obsM) && sameMl(ev, obsB)));
      const ZEvent *realObservation = (realObs.at(i) == -1)
        ? parentTrace.graph.initial() : &(trace.at(realObs.at(i)));
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
