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
(std::vector<ZEvent>&& extension, bool some_to_ann, bool assume_blocked)
  : trace(std::move(extension)), something_to_annotate(some_to_ann),
    has_assume_blocked_thread(assume_blocked), has_error(false)
{
  assert(extension.empty() && !trace.empty());
}


ZExplorer::TraceExtension::TraceExtension
(std::pair<std::vector<ZEvent>&&, bool>&& ext_sometoann)
  : TraceExtension(std::move(ext_sometoann.first),
                   ext_sometoann.second, false) {}


ZExplorer::~ZExplorer()
{
  delete initial;
  initial = nullptr;
}


ZExplorer::ZExplorer(ZBuilderSC& tb)
  : original_TB(&tb)
{
  if (tb.someThreadAssumeBlocked)
    assume_blocked_thread = 1;

  if (tb.somethingToAnnotate.empty())
    executed_traces_full = 1;
  else
    initial = new ZTrace(std::move(tb.prefix), tb.someThreadAssumeBlocked);
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
  std::cout << "Linearization failed:              " << linearization_failed << "\n";
  std::cout << "Linearization succeeded:           " << linearization_succeeded << "\n";
  std::cout << std::setprecision(2) << std::fixed;
  std::cout << "Linearization branching factor:    " << (double)total_children/total_parents << "\n";
  std::cout << "Time spent on copying:             " << time_copy << "\n";
  std::cout << "Time spent on linearization:       " << time_linearization << "\n";
  std::cout << "Time spent on interpreting:        " << time_interpreter << "\n";
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
  #ifndef NDEBUG
    std::cout << "RUNNING DEBUG VERSION\n";
  #endif

  if (initial)
    return explore_rec(*initial);
  return false;
}


bool ZExplorer::explore_rec(ZTrace& ann_trace)
{
  start_err("exploreRec...");
  if (info) {
    llvm::errs() << "-------------------------exploreRec\n\n";
  }
  assert(global_variables_initialized_with_value_zero(ann_trace.trace));

  auto events_to_mutate = ann_trace.events_to_mutate();
  assert(!events_to_mutate.empty());

  // Unset this when any mutation succeeds
  ann_trace.deadlocked = true;

  // Mutate the events
  for (const auto& ev : events_to_mutate) {
    if (is_read(ev)) {
      ann_trace.deadlocked = false;
      bool error = mutate_read(ann_trace, ev);
      if (error) {
        assert(original_TB && original_TB->error_trace);
        end_err("?a");
        return error;
      }
    }
    else {
      assert(is_lock(ev));
      bool error = mutate_lock(ann_trace, ev);
      if (error) {
        assert(original_TB && original_TB->error_trace);
        end_err("?b");
        return error;
      }
    }
    ann_trace.negative.update
      (ev, ann_trace.graph.thread_sizes_minus_one());
  }

  // Done with this recursion node
  if (ann_trace.deadlocked) {
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

bool ZExplorer::mutate_read(const ZTrace& ann_trace, const ZEvent *read)
{
  start_err("mutateRead...");
  assert(is_read(read));
  if (info) {
    llvm::errs() << "Read to mutate:" << *read << "\n";
  }

  std::set<ZAnn> mutation_candidates = ann_trace.mutation_candidates(read);
  if (mutation_candidates.empty()) {
    // All candidates are ruled out by negative annotation
    ++no_mut_choices;
    if (info) {
      llvm::errs() << "All candidates forbidden\n";
    }
    end_err("0a");
    return false;
  }

  for (auto& mutation : mutation_candidates) {
    ++mutations_considered;

    ZAnnotation mutated_annotation(ann_trace.annotation);
    mutated_annotation.add(read->id(), mutation);
    assert(mutated_annotation.size() == ann_trace.annotation.size() + 1);
    auto mutated_po = ann_trace.graph.copy_po();

    if (info) {
      llvm::errs() << "Mutation:\n" << mutation.to_string() << "\n";
    }

    auto init = std::clock();
    ZClosure pre_closure(mutated_annotation, mutated_po);
    pre_closure.rule_one(read, mutation);
    time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

    bool mutation_follows_current_trace = false;
    ZEventID observed_id = (read->observed_trace_id() < 0
      ? ann_trace.graph.initial()->id()
      : ann_trace.trace.at(read->observed_trace_id()).id());
    for (const ZEventID& good : mutation.goodwrites) {
      if (observed_id == good) {
        mutation_follows_current_trace = true;
        break;
      }
    }
    // It is not true that 'read observes value' iff 'read observes a good write'
    // 1) Read could observe a bad write with same value forbidden by negative
    // 2) If 1) happens a lot of times, and despite this we carry forward the
    //    trace all the time, eventually a read could observe a write that is
    //    not even visible for the read in our partial order anymore
    // We should not proceed with this trace in case when only 'read observes value'
    // because we might end up with a full trace and infeasible annotation

    bool error = close_po
      (ann_trace, read, std::move(mutated_annotation),
       std::move(mutated_po), mutation_follows_current_trace);

    if (error) {
      assert(original_TB && original_TB->error_trace);
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

bool ZExplorer::mutate_lock(const ZTrace& ann_trace, const ZEvent *lock)
{
  start_err("mutateLock...");
  assert(is_lock(lock));

  if (info) {
    llvm::errs() << "Lock to mutate: " << *lock << "\n";
  }

  if (!ann_trace.annotation.location_has_some_lock(lock)) {
    // This lock hasn't been touched before
    ann_trace.deadlocked = false;
    if (ann_trace.negative.forbids_initial(lock)) {
      // Negative annotation forbids initial unlock
      end_err("0a");
      if (info) {
        llvm::errs() << "Negative annotation forbids initial unlock\n";
      }
      return false;
    }

    // Trivially realizable
    ++mutations_considered;
    ZAnnotation mutated_annotation(ann_trace.annotation);
    mutated_annotation.set_last_lock(lock);
    auto mutated_po = ann_trace.graph.copy_po();

    if (info) {
      llvm::errs() << "Not locked before\n";
    }

    bool mutation_follows_current_trace =
      (lock->observed_trace_id() == -1);

    auto res = close_po
      (ann_trace, lock, std::move(mutated_annotation),
       std::move(mutated_po), mutation_follows_current_trace);
    end_err("?a");
    return res;
  }

  // The lock has been touched before
  auto last_lock_obs = ann_trace.annotation.last_lock(lock);
  auto last_unlock = ann_trace.graph.unlock_of_this_lock(last_lock_obs);

  if (!last_unlock) {
    // This lock is currently locked
    // Trivially unrealizable
    end_err("0b");
    if (info) {
      llvm::errs() << "Currently locked\n";
    }
    return false;
  }

  // This lock is currently unlocked by lastUnlock
  assert(last_unlock && is_unlock(last_unlock) && same_ml(lock, last_unlock));
  ann_trace.deadlocked = false;

  if (ann_trace.negative.forbids(lock, last_unlock)) {
    // Negative annotation forbids this unlock
    end_err("0c");
    return false;
  }

  // Realizable
  ++mutations_considered;
  ZAnnotation mutated_annotation(ann_trace.annotation);
  mutated_annotation.set_last_lock(lock);
  auto mutated_po = ann_trace.graph.copy_po();

  if (info) {
    llvm::errs() << "Currently unlocked by: " << *last_unlock << "\n";
  }

  auto init = std::clock();
  ZClosure pre_closure(mutated_annotation, mutated_po);
  pre_closure.rule_one_lock(lock, last_unlock);
  time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

  bool mutation_follows_current_trace =
    (lock->observed_trace_id() == last_unlock->trace_id());

  end_err("?b");
  return close_po
    (ann_trace, lock, std::move(mutated_annotation),
     std::move(mutated_po), mutation_follows_current_trace);
}


/* *************************** */
/* CLOSE PO                    */
/* *************************** */

bool ZExplorer::close_po
(const ZTrace& ann_trace, const ZEvent *read_lock,
 ZAnnotation&& mutated_annotation, ZPartialOrder&& mutated_po,
 bool mutation_follows_current_trace)
{
  start_err("closePO...");
  auto init = std::clock();
  ZClosure closure(mutated_annotation, mutated_po);
  bool closed = closure.close
    (is_lock(read_lock) ? nullptr : read_lock);
  double time = (double)(clock() - init)/CLOCKS_PER_SEC;
  time_closure += time;

  if (!closed) {
    ++closure_failed;
    end_err("0");
    return false;
  }

  ++closure_succeeded;
  if (closure.added_edges == 0) {
    closure_no_edge++;
    time_closure_no_edge += time;
  }
  closure_edges += closure.added_edges;
  closure_iter += closure.iterations;

  auto res = extend_and_recur
    (ann_trace, std::move(mutated_annotation),
     std::move(mutated_po), mutation_follows_current_trace);
  end_err("?");
  return res;
}


/* *************************** */
/* EXTEND AND RECUR            */
/* *************************** */

bool ZExplorer::extend_and_recur
(const ZTrace& parent_trace, ZAnnotation&& mutated_annotation,
 ZPartialOrder&& mutated_po, bool mutation_follows_current_trace)
{
  start_err("extendAndRecur...");
  TraceExtension mutated_trace;
  if (mutation_follows_current_trace) {
    mutated_trace = reuse_trace(parent_trace, mutated_annotation);
  }
  else {
    clock_t init = std::clock();
    ZLinearization linearizer(mutated_annotation, mutated_po, parent_trace.trace);
    std::vector<ZEvent> linear;
    linear = linearizer.linearize();
    end_err("finished linearisation1");
    time_linearization += (double)(clock() - init)/CLOCKS_PER_SEC;
    end_err("finished linearisation2");
    //return false;
    // total_parents += linearizer.num_parents;
    // total_children += linearizer.num_children;
    end_err("finished linearisation3");
    if (linear.empty()) {
      ++linearization_failed;
      end_err("0a");
      return false;
    }
    ++linearization_succeeded;
    assert(linearization_respects_annotation(linear, mutated_annotation,
                                             mutated_po, parent_trace));
    mutated_trace = extend_trace(std::move(linear));
    assert(extension_respects_annotation(mutated_trace.trace, mutated_annotation,
                                         mutated_po, parent_trace));
  }

  executed_traces++;
  if (mutated_trace.has_error) {
    assert(!mutation_follows_current_trace);
    end_err("1");
    return true; // Found an error
  }
  if (mutated_trace.has_assume_blocked_thread) {
    // This recursion subtree of the algorithm will only
    // have traces that violate the same assume-condition
    assume_blocked_thread++;
    //return false;
  }
  if (!mutated_trace.something_to_annotate) {
    // Maximal trace
    executed_traces_full++;
    if (info) {
      llvm::errs() << "FULL\n";
    }
    end_err("0b");
    return false;
  }

  clock_t init = std::clock();
  assert(!mutated_trace.empty() && !mutated_po.empty());
  err_msg("creating extension");
  ZTrace mutated_ZTrace
    (parent_trace,
     std::move(mutated_trace.trace),
     std::move(mutated_annotation),
     std::move(mutated_po),
     mutated_trace.has_assume_blocked_thread);
  err_msg("created extension");
  assert(mutated_trace.empty() && mutated_po.empty() &&
         mutated_annotation.empty());
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;

  auto res = explore_rec(mutated_ZTrace);
  end_err("?");
  return res;
}


/* *************************** */
/* REUSE TRACE                 */
/* *************************** */

ZExplorer::TraceExtension
ZExplorer::reuse_trace
(const ZTrace& parent_trace, const ZAnnotation& mutated_annotation)
{
  start_err("reuseTrace..");
  std::vector<ZEvent> tr;
  tr.reserve(parent_trace.trace.size());
  bool something_to_annotate = false;

  for (const ZEvent& ev : parent_trace.trace) {
    if (!something_to_annotate) {
      if (is_read(ev) && (!mutated_annotation.defines(&ev)))
        something_to_annotate = true;
      else if (is_lock(ev)) {
        if (!parent_trace.graph.has_event(&ev))
          something_to_annotate = true;
        else if (ev.event_id() == // last in its thread
                 (parent_trace.graph)(ev.cpid()).size() - 1
                 && !mutated_annotation.is_last_lock(&ev))
          something_to_annotate = true;
      }
    }

    tr.push_back(ZEvent(ev));
  }

  auto res = TraceExtension(std::move(tr), something_to_annotate,
                            parent_trace.assumeblocked);
  end_err("?");
  return res;
}


/* *************************** */
/* EXTEND TRACE                */
/* *************************** */

ZExplorer::TraceExtension
ZExplorer::extend_trace(std::vector<ZEvent>&& tr)
{
  start_err("extendTrace...");
  assert(original_TB);
  clock_t init = std::clock();
  ZBuilderSC TB(*(original_TB->config), original_TB->M, std::move(tr));
  auto trace_extension = TraceExtension(TB.extendGivenTrace());
  time_interpreter += (double)(clock() - init)/CLOCKS_PER_SEC;
  interpreter_used++;

  if (TB.has_error()) {
    // ERROR FOUND
    original_TB->error_trace = TB.get_trace();
    trace_extension.has_error = true;
  }

  if (TB.someThreadAssumeBlocked)
    trace_extension.has_assume_blocked_thread = true;

  end_err("?a");
  return trace_extension;
}


/* *************************** */
/* RESPECTS ANNOTATION         */
/* *************************** */

bool ZExplorer::extension_respects_annotation
(const std::vector<ZEvent>& trace, const ZAnnotation& annotation,
 const ZPartialOrder& mutated_po, const ZTrace& parent_trace) const
{
  assert(!trace.empty());
  for (const ZEvent& evref : trace) {
    const ZEvent *ev = &evref;
    if (is_read(ev)) {
      if (annotation.defines(ev->id())) {
        const ZAnn& ann = annotation.ann(ev->id());
        bool observes_good_write = false;
        const ZEventID observed_id = (ev->observed_trace_id() < 0
          ? ZEventID(CPid(), -1) // initial event id
          : trace.at(ev->observed_trace_id()).id());
        for (const ZEventID& good : ann.goodwrites)
          if (observed_id == good) {
            assert(!observes_good_write);
            observes_good_write = true;
          }
        if (!observes_good_write) {
          parent_trace.dump();
          llvm::errs() << "Closed annotated partial order\n" << mutated_po << "\n";
          llvm::errs() << "Extension\n";
          dump_trace(trace);
          llvm::errs() << "Full annotation that should be respected\n" << annotation << "\n";
          llvm::errs() << "This read         :::  " << ev->to_string(true) << "\n";
          llvm::errs() << "Has annotation    :::  " << ann.to_string() << "\n";
          llvm::errs() << "Actually observed :::  " << observed_id.to_string() << "\n";
          return false;
        }
      }
    }
  }
  return true;
}


bool ZExplorer::linearization_respects_annotation
(const std::vector<ZEvent>& trace, const ZAnnotation& annotation,
 const ZPartialOrder& mutated_po, const ZTrace& parent_trace) const
{
  assert(!trace.empty());
  // Trace index of the last write happening in the given location
  std::unordered_map<SymAddrSize, int> last_write;
  // Trace index of the last lock + whether given ML is currently locked
  std::unordered_map<SymAddrSize, int> last_lock;
  std::unordered_set<SymAddrSize> locked;
  // ThreadID -> ML -> Writes in store queue of thr for ml
  std::unordered_map
    <int, std::unordered_map
     <SymAddrSize, std::list<int>>> storeQueue;
  // Observation in trace: Read -> id in trace of observed event
  std::unordered_map<int, int> observed_trace_id;

  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = &(trace.at(i));

    if (is_write(ev)) {
      last_write[ev->ml()] = i;
    }

    if (is_read(ev)) {
      if (last_write.count(ev->ml()))
        observed_trace_id.emplace(i, last_write.at(ev->ml()));
      else
        observed_trace_id.emplace(i, -1);
    }

    if (is_lock(ev)) {
      assert(!locked.count(ev->ml()));
      locked.insert(ev->ml());
      last_lock[ev->ml()] = i;
    }

    if (is_unlock(ev)) {
      assert(locked.count(ev->ml()));
      locked.erase(ev->ml());
    }
  }

  // Check whether last locks are consistent with annotation
  for (const auto& ml_last : last_lock) {
    const ZEvent *ev = &(trace.at(ml_last.second));
    assert(annotation.is_last_lock(ev));
    if (!annotation.is_last_lock(ev))
      return false;
  }

  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = &(trace.at(i));
    if (is_read(ev)) {
      assert(annotation.defines(ev->id()) &&
             "All reads of the linearization problem input should have desired good-writes defined.");
      const ZAnn& ann = annotation.ann(ev->id());
      bool observes_good_write = false;
      const ZEventID observed_id = (observed_trace_id.at(i) < 0
        ? ZEventID(CPid(), -1) // initial event id
        : trace.at(observed_trace_id.at(i)).id());
      for (const ZEventID& good : ann.goodwrites)
        if (observed_id == good) {
          assert(!observes_good_write);
          observes_good_write = true;
        }
      if (!observes_good_write) {
        parent_trace.dump();
        llvm::errs() << "Closed annotated partial order\n" << mutated_po << "\n";
        llvm::errs() << "Linearization\n";
        dump_trace(trace);
        llvm::errs() << "Full annotation that should be respected\n" << annotation << "\n";
        llvm::errs() << "This read         :::  " << ev->to_string(true) << "\n";
        llvm::errs() << "Has annotation    :::  " << ann.to_string() << "\n";
        llvm::errs() << "Actually observed :::  " << observed_id.to_string() << "\n";
        return false;
      }
    }
  }
  return true;
}


bool ZExplorer::global_variables_initialized_with_value_zero
(const std::vector<ZEvent>& trace) const
{
  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = &(trace.at(i));
    if (is_read(ev)) {
      // If read observes initial, make sure it observes value 0
      if (ev->observed_trace_id() < 0 && ev->value() != 0) {
        dump_trace(trace);
        llvm::errs() << "This read:  " << ev->to_string(true) << "\n";
        llvm::errs() << "Observes the initial event, but not value 0 "
                     << "(it observes value " << ev->value() << ")\n";
        llvm::errs() << "Please check the benchmark, make sure all global "
                     << "variables are initialized with value 0\n";
        return false;
      }
    }
  }
  return true;
}
