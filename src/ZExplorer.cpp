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
(const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& ext_trace,
 const std::shared_ptr<ZGraph>& ext_graph,
 const std::shared_ptr<ZPartialOrder>& ext_po_full,
 const std::shared_ptr<ZPartialOrder>& ext_po_part,
 bool some_to_ann, bool assume_blocked)
  : trace(ext_trace),
    graph(ext_graph),
    po_full(ext_po_full),
    po_part(ext_po_part),
    something_to_annotate(some_to_ann),
    has_assume_blocked_thread(assume_blocked), has_error(false)
{
  assert(!empty());
  assert(ext_trace && !ext_trace->empty());
  assert(ext_graph && !ext_graph->empty());
  assert(ext_po_full && !ext_po_full->empty());
  assert(ext_po_part && !ext_po_part->empty());
}


bool ZExplorer::TraceExtension::empty() const
{
  return ((!trace || trace->empty()) &&
          (!graph || graph->empty()) &&
          (!po_full || po_full->empty()) &&
          (!po_part || po_part->empty()));
}


ZExplorer::ZExplorer(ZBuilderSC& tb)
  : original_TB(&tb)
{
  if (tb.someThreadAssumeBlocked)
    assume_blocked_thread = 1;

  if (tb.somethingToAnnotate.empty())
    executed_traces_full = 1;
  else
    initial = std::unique_ptr<ZTrace>(new ZTrace
    (tb.prefix, tb.graph, tb.po_full, tb.po_part, tb.someThreadAssumeBlocked));
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
  double avg1_branch = (total_parents==0) ? 1.0 : ((double)total_children/total_parents);
  std::cout << "Linearization branching Avg1:      " << avg1_branch << "\n";
  std::cout << "Linearization branching Avg2:      " << avg2_branch << "\n";
  std::cout << "Linearization branching Max:       " << max_branch << "\n";
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
  start_err("explore_rec...");
  assert(global_variables_initialized_with_value_zero(ann_trace.trace()));

  auto events_to_mutate = ann_trace.events_to_mutate();
  assert(!events_to_mutate.empty());

  // Unset this when any mutation succeeds
  ann_trace.deadlocked = true;

  auto negative_update = ann_trace.po_part().thread_sizes_minus_one();

  // Collect mutations
  for (const auto& ev : events_to_mutate) {
    if (is_read(ev)) {
      ann_trace.deadlocked = false;
      bool error = mutate_read(ann_trace, ev);
      if (error) {
        assert(original_TB && original_TB->error_trace);
        end_err("1a");
        return error;
      }
    }
    else {
      assert(is_lock(ev));
      bool error = mutate_lock(ann_trace, ev);
      if (error) {
        assert(original_TB && original_TB->error_trace);
        end_err("1b");
        return error;
      }
    }
  }

  if (ann_trace.deadlocked) {
    // Not full trace ending in a deadlock
    // We count it as full
    assert(ann_trace.children_lock.empty());
    assert(ann_trace.children_read.empty());
    ++executed_traces_full;
    ++executed_traces_full_deadlock;
  }

  // Recursive calls - locks
  for (const auto& lock : events_to_mutate) {
    if (!is_lock(lock))
      continue;
    if (!ann_trace.children_lock.count(lock))
      continue;
    ZTrace& mutated_trace = *ann_trace.children_lock.at(lock).get();
    mutated_trace.set_negative(ann_trace.negative());
    bool error = explore_rec(mutated_trace);
    if (error) {
      assert(original_TB && original_TB->error_trace);
      end_err("1c");
      return error;
    }
    ann_trace.negative().update(lock, negative_update);
  }

  // Recursive calls - reads
  for (const auto& read : events_to_mutate) {
    if (!is_read(read))
      continue;
    if (!ann_trace.children_read.count(read))
      continue;
    for (const auto& value_mutatedtrace : ann_trace.children_read.at(read)) {
      ZTrace& mutated_trace = *value_mutatedtrace.second.get();
      mutated_trace.set_negative(ann_trace.negative());
      bool error = explore_rec(mutated_trace);
      if (error) {
        assert(original_TB && original_TB->error_trace);
        end_err("1d");
        return error;
      }
    }
    ann_trace.negative().update(read, negative_update);
  }

  // Done with this recursion node
  end_err("0");
  return false;
}


/* *************************** */
/* MUTATE READ                 */
/* *************************** */

bool ZExplorer::mutate_read(ZTrace& ann_trace, const ZEvent *read)
{
  start_err("mutate_read...");
  assert(is_read(read));

  std::set<ZAnn> mutation_candidates = ann_trace.mutation_candidates(read);
  if (mutation_candidates.empty()) {
    // All candidates are ruled out by negative annotation
    ++no_mut_choices;
    end_err("0a");
    return false;
  }

  for (auto& mutation : mutation_candidates) {
    err_msg(mutation.to_string());
    ++mutations_considered;

    ZAnnotation mutated_annotation(ann_trace.annotation());
    mutated_annotation.add(read->id(), mutation);
    assert(mutated_annotation.size() == ann_trace.annotation().size() + 1);
    ZPartialOrder mutated_po(ann_trace.po_part());

    auto init = std::clock();
    ZClosure pre_closure(mutated_annotation, mutated_po);
    pre_closure.rule_one(read, mutation);
    time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

    bool mutation_follows_current_trace = false;
    assert(read->observed_trace_id() < (int) ann_trace.trace().size());
    ZEventID observed_id = (read->observed_trace_id() < 0
      ? ann_trace.graph().initial()->id()
      : ann_trace.trace().at(read->observed_trace_id())->id());
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

    auto founderror_mutatedtrace = close_po
      (ann_trace, read, std::move(mutated_annotation),
       std::move(mutated_po), mutation_follows_current_trace);

    if (founderror_mutatedtrace.first) {
      assert(original_TB && original_TB->error_trace);
      end_err("1");
      return true;
    }
    if (!founderror_mutatedtrace.second)
      continue;
    if (!ann_trace.children_read.count(read)) {
      ann_trace.children_read.emplace(
        read, std::map<int, std::unique_ptr<ZTrace>>());
    }
    ann_trace.children_read[read].emplace(
      mutation.value, std::move(founderror_mutatedtrace.second));
  }

  end_err("0b");
  return false;
}


/* *************************** */
/* MUTATE LOCK                 */
/* *************************** */

bool ZExplorer::mutate_lock(ZTrace& ann_trace, const ZEvent *lock)
{
  start_err("mutate_lock...");
  assert(is_lock(lock));

  if (!ann_trace.annotation().location_has_some_lock(lock)) {
    // This lock hasn't been touched before
    ann_trace.deadlocked = false;
    if (ann_trace.negative().forbids_initial(lock)) {
      // Negative annotation forbids initial unlock
      end_err("0a");
      return false;
    }

    // Trivially realizable
    ++mutations_considered;
    ZAnnotation mutated_annotation(ann_trace.annotation());
    mutated_annotation.set_last_lock(lock);
    ZPartialOrder mutated_po(ann_trace.po_part());

    bool mutation_follows_current_trace =
      (lock->observed_trace_id() == -1);

    auto founderror_mutatedtrace = close_po
      (ann_trace, lock, std::move(mutated_annotation),
       std::move(mutated_po), mutation_follows_current_trace);
    if (founderror_mutatedtrace.first) {
      end_err("1a");
      return true;
    }
    assert(!ann_trace.children_lock.count(lock));
    if (founderror_mutatedtrace.second) {
      ann_trace.children_lock.emplace
      (lock, std::move(founderror_mutatedtrace.second));
    }
    end_err("?a");
    return false;
  }

  // The lock has been touched before
  auto last_lock_obs = ann_trace.annotation().last_lock(lock);
  auto last_unlock = ann_trace.graph().unlock_of_this_lock(last_lock_obs);

  if (!last_unlock || !ann_trace.po_part().spans_event(last_unlock)) {
    // This lock is currently locked
    // Trivially unrealizable
    end_err("0b");
    return false;
  }

  // This lock is currently unlocked by lastUnlock
  assert(last_unlock && is_unlock(last_unlock) && same_ml(lock, last_unlock));
  assert(ann_trace.po_part().spans_event(last_unlock));
  ann_trace.deadlocked = false;

  if (ann_trace.negative().forbids(lock, last_unlock)) {
    // Negative annotation forbids this unlock
    end_err("0c");
    return false;
  }

  // Realizable, currently unlocked by last_unlock
  ++mutations_considered;
  ZAnnotation mutated_annotation(ann_trace.annotation());
  mutated_annotation.set_last_lock(lock);
  ZPartialOrder mutated_po(ann_trace.po_part());

  auto init = std::clock();
  ZClosure pre_closure(mutated_annotation, mutated_po);
  pre_closure.rule_one_lock(lock, last_unlock);
  time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

  bool mutation_follows_current_trace =
    (lock->observed_trace_id() == last_unlock->trace_id());

  auto founderror_mutatedtrace = close_po
    (ann_trace, lock, std::move(mutated_annotation),
     std::move(mutated_po), mutation_follows_current_trace);
  if (founderror_mutatedtrace.first) {
    end_err("1b");
    return true;
  }
  assert(!ann_trace.children_lock.count(lock));
  if (founderror_mutatedtrace.second) {
    ann_trace.children_lock.emplace
    (lock, std::move(founderror_mutatedtrace.second));
  }
  end_err("?b");
  return false;
}


/* *************************** */
/* CLOSE PO                    */
/* *************************** */

std::pair<bool, std::unique_ptr<ZTrace>> ZExplorer::close_po
(const ZTrace& ann_trace, const ZEvent *read_lock,
 ZAnnotation&& mutated_annotation, ZPartialOrder&& mutated_po,
 bool mutation_follows_current_trace)
{
  start_err("close_po...");
  auto init = std::clock();
  ZClosure closure(mutated_annotation, mutated_po);
  bool closed = closure.close
    (is_lock(read_lock) ? nullptr : read_lock);
  double time = (double)(clock() - init)/CLOCKS_PER_SEC;
  time_closure += time;

  if (!closed) {
    ++closure_failed;
    end_err("0");
    return {false, std::unique_ptr<ZTrace>(nullptr)};
  }

  ++closure_succeeded;
  if (closure.added_edges == 0) {
    closure_no_edge++;
    time_closure_no_edge += time;
  }
  closure_edges += closure.added_edges;
  closure_iter += closure.iterations;

  auto res = realize_mutation
    (ann_trace, read_lock, std::move(mutated_annotation),
     std::move(mutated_po), mutation_follows_current_trace);
  end_err("?");
  return res;
}


/* *************************** */
/* EXTEND AND RECUR            */
/* *************************** */

std::pair<bool, std::unique_ptr<ZTrace>> ZExplorer::realize_mutation
(const ZTrace& parent_trace, const ZEvent *read_lock,
 ZAnnotation&& mutated_annotation, ZPartialOrder&& mutated_po,
 bool mutation_follows_current_trace)
{
  start_err("realize_mutation...");
  TraceExtension mutated_trace;
  if (mutation_follows_current_trace) {
    mutated_trace = reuse_trace
    (parent_trace, read_lock, mutated_annotation, std::move(mutated_po));
  }
  else {
    clock_t init = std::clock();
    ZLinearization linearizer(mutated_annotation, mutated_po);
    std::vector<ZEvent> linear = linearizer.linearize();
    time_linearization += (double)(clock() - init)/CLOCKS_PER_SEC;
    ++total_lin;
    total_parents += linearizer.num_parents;
    total_children += linearizer.num_children;
    double cur_br = (linearizer.num_parents==0) ? 1.0 :
      ((double)linearizer.num_children/linearizer.num_parents);
    avg2_branch += ((double)(cur_br - avg2_branch)/total_lin);
    if (cur_br > max_branch) {
      max_branch = cur_br;
    }
    err_msg("finished linearisation");
    if (linear.empty()) {
      ++linearization_failed;
      end_err("0a-linfailed");
      return {false, std::unique_ptr<ZTrace>(nullptr)};
    }
    ++linearization_succeeded;
    assert(total_lin == linearization_succeeded + linearization_failed);
    assert(linearization_respects_annotation(linear, mutated_annotation,
                                             mutated_po, parent_trace));
    assert(!mutated_po.empty() && !linear.empty());
    mutated_trace = extend_trace(std::move(linear), std::move(mutated_po));
    assert(mutated_po.empty() && linear.empty());
    assert(extension_respects_annotation(*mutated_trace.trace, mutated_annotation,
                                         mutated_po, parent_trace));
  }

  executed_traces++;
  if (mutated_trace.has_error) {
    assert(!mutation_follows_current_trace);
    end_err("1-error");
    return {true, std::unique_ptr<ZTrace>(nullptr)}; // Found an error
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
    end_err("0b-full");
    return {false, std::unique_ptr<ZTrace>(nullptr)};
  }

  clock_t init = std::clock();
  assert(mutated_po.empty());
  err_msg("creating extension");
  std::unique_ptr<ZTrace> mutated_ZTrace(new ZTrace(
    parent_trace,
    mutated_trace.trace,
    std::move(mutated_annotation),
    mutated_trace.graph,
    mutated_trace.po_full,
    mutated_trace.po_part,
    mutated_trace.has_assume_blocked_thread));
  err_msg("created extension");
  assert(mutated_annotation.empty() && mutated_po.empty());
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;

  end_err("?-realized");
  return {false, std::move(mutated_ZTrace)};
}


/* *************************** */
/* REUSE TRACE                 */
/* *************************** */

ZExplorer::TraceExtension
ZExplorer::reuse_trace
(const ZTrace& parent_trace, const ZEvent *read_lock,
 const ZAnnotation& mutated_annotation, ZPartialOrder&& mutated_po)
{
  start_err("reuse_trace...");
  bool something_to_annotate = false;
  for (int i = parent_trace.trace().size() - 1; i >= 0; --i) {
    const ZEvent& ev = *(parent_trace.trace().at(i));
    assert(parent_trace.graph().has_event(&ev));
    assert(!something_to_annotate);
    if (is_read(ev) && (!mutated_annotation.defines(&ev))) {
      something_to_annotate = true;
      break;
    }
    if (is_lock(ev)) {
      if (!parent_trace.po_part().spans_event(&ev)) {
        something_to_annotate = true;
        break;
      }
      else if (ev.event_id() == // last in its thread
               parent_trace.po_part().thread_size(ev.cpid()) - 1
               && !mutated_annotation.is_last_lock(&ev)) {
        something_to_annotate = true;
        break;
      }
    }
  }

  clock_t init = std::clock();
  assert(!mutated_po.empty());
  std::shared_ptr<ZPartialOrder> mutated_po_ptr(new ZPartialOrder
  (std::move(mutated_po), parent_trace.graph()));
  assert(mutated_po.empty());
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;

  if (something_to_annotate) {
    mutated_po_ptr->extend(read_lock, mutated_annotation);
  }

  init = std::clock();
  auto res = TraceExtension
  (parent_trace.trace_ptr(),
   parent_trace.graph_ptr(),
   parent_trace.po_full_ptr(),
   mutated_po_ptr,
   something_to_annotate,
   parent_trace.assumeblocked);
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;

  end_err("?");
  return res;
}


/* *************************** */
/* EXTEND TRACE                */
/* *************************** */

ZExplorer::TraceExtension
ZExplorer::extend_trace
(std::vector<ZEvent>&& tr, ZPartialOrder&& mutated_po)
{
  start_err("extend_trace...");
  assert(original_TB);
  assert(!tr.empty() && !mutated_po.empty());
  clock_t init = std::clock();
  ZBuilderSC TB(*(original_TB->config), original_TB->M,
                std::move(tr), std::move(mutated_po));
  assert(tr.empty() && mutated_po.empty());
  TB.extendGivenTrace();
  time_interpreter += (double)(clock() - init)/CLOCKS_PER_SEC;
  interpreter_used++;
  assert(TB.prefix && !TB.prefix->empty());
  assert(TB.graph && !TB.graph->empty());
  assert(TB.po_full && !TB.po_full->empty());
  assert(TB.po_part && !TB.po_part->empty());

  init = std::clock();
  TraceExtension trace_extension(
    TB.prefix, TB.graph, TB.po_full, TB.po_part,
    !TB.somethingToAnnotate.empty(), TB.someThreadAssumeBlocked);
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;

  if (TB.has_error()) {
    // ERROR FOUND
    original_TB->error_trace = TB.get_trace();
    trace_extension.has_error = true;
  }

  end_err("?a");
  return trace_extension;
}


/* *************************** */
/* RESPECTS ANNOTATION         */
/* *************************** */

bool ZExplorer::extension_respects_annotation
(const std::vector<std::unique_ptr<ZEvent>>& trace,
 const ZAnnotation& annotation, const ZPartialOrder& mutated_po,
 const ZTrace& parent_trace) const
{
  assert(!trace.empty());
  for (const std::unique_ptr<ZEvent>& evref : trace) {
    const ZEvent *ev = evref.get();
    if (is_read(ev)) {
      if (annotation.defines(ev->id())) {
        const ZAnn& ann = annotation.ann(ev->id());
        bool observes_good_write = false;
        assert(ev->observed_trace_id() < (int) trace.size());
        const ZEventID observed_id = (ev->observed_trace_id() < 0
          ? ZEventID(CPid(), -1) // initial event id
          : trace.at(ev->observed_trace_id())->id());
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
(const std::vector<std::unique_ptr<ZEvent>>& trace) const
{
  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = trace.at(i).get();
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
