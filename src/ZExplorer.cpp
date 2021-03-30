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
  std::cout << "Early stopping failed:             " << early_failed << "\n";
  std::cout << "Early stopping succeeded:          " << early_succeeded << "\n";
  std::cout << "Reads with no mutation choices:    " << no_mut_choices << "\n";
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
  std::cout << "Time linearization noaux:          " << time_lin_noaux << "\n";
  std::cout << "Time linearization noclos:         " << time_lin_nocl << "\n";
  std::cout << "Time linearization noclos noaux:   " << time_lin_noclnoaux << "\n";
  std::cout << "Time spent on interpreting:        " << time_interpreter << "\n";
  std::cout << "Time spent on closure:             " << time_closure << "\n";
  std::cout << "Time spent on closure-succ-noedge: " << time_closure_no_edge << "\n";
  std::cout << "Time spent on early stopping:      " << time_early << "\n";
  std::cout << "\n" << std::scientific;

  // Change to false to test if assertions are on
  // To disable assertions (i.e. build as Release),
  // in src/Makefile add in CXXFLAGS this: -DNDEBUG
  assert(true && "RUN ON RELEASE");
}

/* *************************** */
/* ORDER TO MUTATE             */
/* *************************** */

std::vector<const ZEvent *> ZExplorer::get_order_to_mutate
(const ZTrace& ann_trace,
 const std::list<const ZEvent *>& locks_to_mutate,
 const std::list<const ZEvent *>& reads_to_mutate_noneg,
 const std::list<const ZEvent *>& reads_to_mutate_withneg) const
{
  std::vector<const ZEvent *> result;
  result.reserve(locks_to_mutate.size() + reads_to_mutate_noneg.size() +
                 reads_to_mutate_withneg.size());

  for (const auto& lock : locks_to_mutate) {
    assert(is_lock(lock));
    result.push_back(lock);
  }

  get_order_to_mutate_help(ann_trace, reads_to_mutate_noneg, result);

  for (const auto& read : reads_to_mutate_withneg) {
    assert(is_read(read));
    result.push_back(read);
  }

  assert(result.size() == locks_to_mutate.size() +
         reads_to_mutate_noneg.size() + reads_to_mutate_withneg.size());
  return result;
}


void ZExplorer::get_order_to_mutate_help
(const ZTrace& ann_trace,
 const std::list<const ZEvent *>& reads_to_mutate_noneg,
 std::vector<const ZEvent *>& order) const
{
  const std::unordered_map
  <const ZEvent *, std::unordered_set<SymAddrSize>>& unc = (
    ann_trace.graph().cache().read_uncovers_mls);
  std::list<const ZEvent *> uncovers_nothing;
  std::set<const ZEvent *, ZEventPtrComp> greedy;
  std::unordered_map<SymAddrSize, int> how_many_reads_of_ml;
  for (const auto& read : reads_to_mutate_noneg) {
    assert(is_read(read));
    assert(unc.count(read));
    if (unc.at(read).empty()) {
      uncovers_nothing.push_front(read);
      continue;
    }
    assert(!unc.at(read).empty());
    greedy.insert(read);
    if (!how_many_reads_of_ml.count(read->ml())) {
      how_many_reads_of_ml.emplace(read->ml(), 1);
    } else {
      how_many_reads_of_ml[read->ml()]++;
    }
  }
#ifndef NDEBUG
  int total = 0;
  for (const auto& ml_howmany : how_many_reads_of_ml) {
    total += ml_howmany.second;
  }
  assert(total == greedy.size());
#endif

  while (!greedy.empty()) {
    // Select the read that uncovers the most mls for the others.
    // That means: for each ml with >=1 write to that ml after read,
    // bump the score by how many other reads of greedy read from that ml.
    const ZEvent * best = nullptr;
    int bestscore = -1;
    //
    for (const auto& read : greedy) {
      int score = (unc.at(read).count(read->ml())) ? -1 : 0;
      for (const SymAddrSize& write_ml_after : unc.at(read)) {
        if (how_many_reads_of_ml.count(write_ml_after)) {
          score += how_many_reads_of_ml[write_ml_after];
        }
      }
      assert(score >= 0);
      if (score > bestscore) {
        best = read;
        bestscore = score;
      }
    }
    //
    assert(is_read(best));
    order.push_back(best);
    how_many_reads_of_ml[best->ml()]--;
    greedy.erase(best);
  }
#ifndef NDEBUG
  for (const auto& ml_howmany : how_many_reads_of_ml) {
    assert(ml_howmany.second == 0);
  }
#endif

  for (const auto& read : uncovers_nothing) {
    assert(is_read(read));
    order.push_back(read);
  }
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

  std::list<const ZEvent *> events_to_mutate = ann_trace.events_to_mutate();
  assert(!events_to_mutate.empty());

  // Unset this when any mutation succeeds
  ann_trace.deadlocked = true;

  auto negative_update = ann_trace.po_part().thread_sizes_minus_one();

  // Collect mutations
  std::map<const ZEvent *, std::set<ZAnn>> read_mutations;
  std::map<const ZEvent *, const ZEvent *> lock_mutations;
  std::list<const ZEvent *> locks_to_mutate;
  std::list<const ZEvent *> reads_to_mutate_noneg;
  std::list<const ZEvent *> reads_to_mutate_withneg;
  for (const auto& ev : events_to_mutate) {
    assert(!read_mutations.count(ev) && !lock_mutations.count(ev));
    if (is_read(ev)) {
      // Collect read mutations
      ann_trace.deadlocked = false;
      read_mutations.emplace(ev, ann_trace.mutation_candidates(ev));
      if (!read_mutations[ev].empty()) {
        // This read has mutation(s) to consider
        if (ann_trace.negative().forbids_initial(ev)) {
          reads_to_mutate_withneg.push_back(ev);
        } else {
          reads_to_mutate_noneg.push_back(ev);
        }
      }
      continue;
    }
    // Collect lock mutation resp. which thread blocks it
    assert(is_lock(ev));
    const ZEvent * mut = collect_lock_mutation(ann_trace, ev);
    if (!mut) {
      // We cannot mutate because negative annotation forbids
      continue;
    }
    // This lock has mutation to consider
    assert(is_initial(mut) || is_unlock(mut));
    lock_mutations.emplace(ev, mut);
    locks_to_mutate.push_back(ev);
  }

  if (info) {
    ann_trace.annotation().dump();
  }

  if (ann_trace.deadlocked) {
    // Not full trace ending in a deadlock
    // We count it as full
    ++executed_traces_full;
    ++executed_traces_full_deadlock;
    end_err("0-deadlock");
    return false;
  }

  if (reads_to_mutate_noneg.empty() && reads_to_mutate_withneg.empty() &&
      locks_to_mutate.empty()) {
    // All mutations forbidden by negative annotation
    end_err("0-allforbidden");
    return false;
  }

  // Early stopping, disabled until its TODO is solved
  auto init = std::clock();
  bool stop_early = early_stopping(ann_trace, read_mutations, lock_mutations);
  time_early += (double)(clock() - init)/CLOCKS_PER_SEC;
  if (stop_early) {
    // This exploration provably leads to no new behaviour, stop early
    end_err("0-early");
    return false;
  }

  std::vector<const ZEvent *> order_to_mutate = get_order_to_mutate(
    ann_trace, locks_to_mutate, reads_to_mutate_noneg, reads_to_mutate_withneg);

  // Recursive calls, check for backtrack signal
  for (int idx = 0; idx < order_to_mutate.size(); ++idx) {
    const ZEvent * read_lock = order_to_mutate[idx];
    assert(is_read(read_lock) || is_lock(read_lock));
    // Initialize backtrack signal
    if (is_read(read_lock)) {
      ann_trace.backtrack = (idx >= locks_to_mutate.size() + reads_to_mutate_noneg.size());
      assert(ann_trace.backtrack || !ann_trace.negative().forbids_initial(read_lock));
      assert(!ann_trace.backtrack || ann_trace.negative().forbids_initial(read_lock));
    } else {
      assert(is_lock(read_lock));
      ann_trace.backtrack = false;
      for (int l_idx = idx + 1; l_idx < order_to_mutate.size(); ++l_idx) {
        if (is_read(order_to_mutate[l_idx])) {
          // Remaining are reads
          #ifndef NDEBUG
          for (int rem = l_idx + 1; rem < order_to_mutate.size(); ++rem) {
            assert(is_read(order_to_mutate[rem]));
          }
          #endif
          break;
        }
        assert(is_lock(order_to_mutate[l_idx]));
        if (same_ml(read_lock, order_to_mutate[l_idx])) {
          // His unlock will be a new source for read_lock
          ann_trace.backtrack = true;
          break;
        }
      }
    }
    // Recursive calls
    if (is_read(read_lock)) {
      assert(read_mutations.count(read_lock) && !read_mutations[read_lock].empty());
      for (const ZAnn& mutation : read_mutations[read_lock]) {
        bool error = mutate_read(ann_trace, read_lock, mutation);
        if (error) {
          assert(original_TB && original_TB->error_trace);
          end_err("1-read");
          return error;
        }
      }
      ann_trace.negative().update(read_lock, negative_update);
    } else {
      assert(lock_mutations.count(read_lock) && lock_mutations[read_lock]);
      bool error = mutate_lock(ann_trace, read_lock, lock_mutations[read_lock]);
      if (error) {
        assert(original_TB && original_TB->error_trace);
        end_err("1-lock");
        return error;
      }
      ann_trace.negative().update(read_lock, negative_update);
    }
    // Done if no backtrack signal
    if (!ann_trace.backtrack) {
      break;
    }
  }

  // Done with this recursion node
  end_err("0");
  return false;
}


/* *************************** */
/* MUTATE READ                 */
/* *************************** */

bool ZExplorer::mutate_read
(ZTrace& ann_trace, const ZEvent *read, const ZAnn& mutation)
{
  start_err(std::string("mutate_read...") + mutation.to_string());
  assert(is_read(read));
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
  assert(read->observed_trace_id() >= -1);
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

  bool error = close_po
    (ann_trace, read, std::move(mutated_annotation),
     std::move(mutated_po), mutation_follows_current_trace);

  if (error) {
    assert(original_TB && original_TB->error_trace);
    end_err("1");
    return true;
  }

  end_err("0");
  return false;
}


/* *************************** */
/* COLLECT LOCK MUTATION       */
/* *************************** */

const ZEvent * ZExplorer::collect_lock_mutation
(const ZTrace& ann_trace, const ZEvent *lock)
{
  start_err("collect_lock_mutation...");
  assert(is_lock(lock));

  if (!ann_trace.annotation().location_has_some_lock(lock)) {
    // This lock hasn't been touched before
    ann_trace.deadlocked = false;
    if (ann_trace.negative().forbids_initial(lock)) {
      // Negative annotation forbids initial unlock
      end_err("0a");
      return nullptr;
    }

    // Trivially realizable
    end_err("?a");
    return ann_trace.graph().initial();
  }

  // The lock has been touched before
  const ZEventID& last_lock_obs = ann_trace.annotation().last_lock(lock);
  const ZEvent *last_unlock = ann_trace.graph().unlock_of_this_lock(last_lock_obs);

  if (!last_unlock || !ann_trace.po_part().spans_event(last_unlock)) {
    // This lock is currently locked
    // Trivially unrealizable
    const ZEvent * holds_lock = ann_trace.graph().event(last_lock_obs);
    assert(is_lock(holds_lock) && ann_trace.po_part().spans_event(holds_lock));
    end_err("0b");
    return nullptr;
  }

  // This lock is currently unlocked by last_unlock
  assert(last_unlock && is_unlock(last_unlock) && same_ml(lock, last_unlock));
  assert(ann_trace.po_part().spans_event(last_unlock));
  ann_trace.deadlocked = false;

  if (ann_trace.negative().forbids(lock, last_unlock)) {
    // Negative annotation forbids this unlock
    end_err("0c");
    return nullptr;
  }

  // Realizable, currently unlocked by last_unlock
  end_err("?b");
  return last_unlock;
}


/* *************************** */
/* MUTATE LOCK                 */
/* *************************** */

bool ZExplorer::mutate_lock
(ZTrace& ann_trace, const ZEvent *lock, const ZEvent *unlock)
{
  start_err("mutate_lock...");
  assert(is_lock(lock));
  assert(is_initial(unlock) || is_unlock(unlock));
  ++mutations_considered;

  if (is_initial(unlock)) {
    // Trivially realizable initial acquire of the lock
    ZAnnotation mutated_annotation(ann_trace.annotation());
    mutated_annotation.set_last_lock(lock);
    ZPartialOrder mutated_po(ann_trace.po_part());

    bool mutation_follows_current_trace =
      (lock->observed_trace_id() == -1);

    bool error = close_po
      (ann_trace, lock, std::move(mutated_annotation),
       std::move(mutated_po), mutation_follows_current_trace);
    if (error) {
      end_err("1a");
      return true;
    }
    end_err("0a");
    return false;
  }

  assert(is_unlock(unlock));
  // Realizable acquire after releasing by unlock
  ZAnnotation mutated_annotation(ann_trace.annotation());
  mutated_annotation.set_last_lock(lock);
  ZPartialOrder mutated_po(ann_trace.po_part());

  auto init = std::clock();
  ZClosure pre_closure(mutated_annotation, mutated_po);
  pre_closure.rule_one_lock(lock, unlock);
  time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;

  bool mutation_follows_current_trace =
    (lock->observed_trace_id() == unlock->trace_id());

  bool error = close_po
    (ann_trace, lock, std::move(mutated_annotation),
     std::move(mutated_po), mutation_follows_current_trace);
  if (error) {
    end_err("1b");
    return true;
  }
  end_err("0b");
  return false;
}


/* *************************** */
/* CLOSE PO                    */
/* *************************** */

bool ZExplorer::close_po
(ZTrace& ann_trace, const ZEvent *read_lock,
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
    return false;
  }

  ++closure_succeeded;
  if (closure.added_edges == 0) {
    closure_no_edge++;
    time_closure_no_edge += time;
  }
  closure_edges += closure.added_edges;
  closure_iter += closure.iterations;

  bool res = realize_mutation
    (ann_trace, read_lock, std::move(mutated_annotation),
     std::move(mutated_po), mutation_follows_current_trace);
  end_err("?");
  return res;
}


/* *************************** */
/* REALIZE MUTATION            */
/* *************************** */

bool ZExplorer::realize_mutation
(ZTrace& parent_trace, const ZEvent *read_lock,
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
    // Prepare thread order + mutex edges
    ZPartialOrder thrord_mutexedges(mutated_po, parent_trace.po_full());
    // Linearization + closure + auxiliary
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
      end_err("0-linfailed");
      return false;
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
    if (info) {
      mutated_annotation.dump();
    }
    executed_traces_full++;
    end_err("0-full");
    return false;
  }

  clock_t init = std::clock();
  assert(mutated_po.empty());
  err_msg("creating extension");
  ZTrace mutated_ZTrace(
    parent_trace,
    mutated_trace.trace,
    std::move(mutated_annotation),
    mutated_trace.graph,
    mutated_trace.po_full,
    mutated_trace.po_part,
    mutated_trace.has_assume_blocked_thread);
  err_msg("created extension");
  assert(mutated_annotation.empty() && mutated_po.empty());
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;

  if (!parent_trace.backtrack) {
    // Record as ancestor to look for backtrack signal
    if (!ancestors.count(read_lock->ml())) {
      ancestors.emplace(read_lock->ml(), std::map<ZEventID, ZTrace *>());
    }
    assert(!ancestors[read_lock->ml()].count(read_lock->id()));
    ancestors[read_lock->ml()].emplace(read_lock->id(), &parent_trace);
  }
  // Recursive call
  bool error = explore_rec(mutated_ZTrace);
  if (ancestors.count(read_lock->ml()) &&
      ancestors[read_lock->ml()].count(read_lock->id())) {
    // Unrecord as ancestor
    assert(ancestors[read_lock->ml()][read_lock->id()] == &parent_trace);
    ancestors[read_lock->ml()].erase(read_lock->id());
    if (ancestors[read_lock->ml()].empty()) {
      ancestors.erase(read_lock->ml());
    }
  }

  end_err("?-realized");
  return error;
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
  assert((is_read(read_lock) && mutated_annotation.defines(read_lock)) ||
         (is_lock(read_lock) && mutated_annotation.is_last_lock(read_lock)));

  clock_t init = std::clock();
  bool something_to_annotate = (
    parent_trace.po_part().something_to_annotate(mutated_annotation));

  assert(!mutated_po.empty());
  std::shared_ptr<ZPartialOrder> mutated_po_ptr(new ZPartialOrder
  (std::move(mutated_po), parent_trace.graph()));
  assert(mutated_po.empty());

  if (something_to_annotate) {
    mutated_po_ptr->extend(read_lock, mutated_annotation,
                           *this, parent_trace.po_full());
  } else {
    mutated_po_ptr->process_remaining_events_for_backtrack_points(
                    *this, parent_trace.po_full());
  }

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
                std::move(tr), std::move(mutated_po), this);
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
/* BACKTRACK POINTS            */
/* *************************** */


bool ZExplorer::try_add_backtrack_point
(const ZPartialOrder& po_full, const ZEvent * new_src,
 ZTrace * anc_trace, const ZEventID& ancid) const
{
  assert(anc_trace);
  assert(!anc_trace->backtrack);
  // Disregard if same thread
  if (ancid.cpid() == new_src->cpid()) {
    assert(ancid.event_id() < new_src->event_id());
    return false;
  }
  // Disregard if ancev --thread_order--> new_src
  const ZEvent * ancev = po_full.graph.event(ancid);
  assert(same_ml(ancev, new_src));
  assert((is_read(ancev) && is_write(new_src)) ||
         (is_lock(ancev) && is_lock(new_src)) ||
         (is_read(ancev) && new_src->is_read_of_cas()));
  assert(po_full.spans_event(ancev));
  assert(!po_full.has_edge(new_src, ancev));
  if (po_full.has_edge(ancev, new_src)) {
    return false;
  }
  // TODO: Disregard if given the closure edges so far in anc_trace->po_part,
  // new_src will surely not be observable (e.g. ancev -> new_src will happen)
  // Add backtrack point
  anc_trace->backtrack = true;
  return true;
}


void ZExplorer::process_backtrack_points
(const ZPartialOrder& po_full, const ZEvent * new_src) const
{
  assert(is_write(new_src) || is_lock(new_src) || new_src->is_read_of_cas());
  assert(po_full.graph.has_event(new_src));
  assert(po_full.spans_event(new_src));
  assert(ancestors.count(new_src->ml()));
  for (auto it = ancestors[new_src->ml()].begin();
            it != ancestors[new_src->ml()].end(); ) {
    const ZEventID& ancid = it->first;
    assert(po_full.graph.has_event(ancid));
    ZTrace * anc_trace = it->second;
    assert(anc_trace);
    bool added = try_add_backtrack_point(
      po_full, new_src, anc_trace, ancid);
    if (added) {
      assert(anc_trace->backtrack);
      it = ancestors[new_src->ml()].erase(it);
    } else {
      assert(!anc_trace->backtrack);
      ++it;
    }
  }
  if (ancestors[new_src->ml()].empty()) {
    ancestors.erase(new_src->ml());
  }
}


/* *************************** */
/* EARLY STOPPING              */
/* *************************** */

bool ZExplorer::early_stopping
(const ZTrace& ann_trace,
 const std::map<const ZEvent *, std::set<ZAnn>>& read_mutations,
 const std::map<const ZEvent *, const ZEvent *>& lock_mutations)
{
  // If there is a read with multiple values to see, no early stopping
  bool read_with_more = false;
  std::set<const ZEvent *> maybe_reads;
  for (const auto& r_mut : read_mutations) {
    assert(is_read(r_mut.first));
    if (r_mut.second.empty()) {
      // All mutation candidates are ruled out by negative annotation
      ++no_mut_choices;
      maybe_reads.insert(r_mut.first);
    } else if (r_mut.second.size() > 1) {
      // Read can see different values, no early stopping
      read_with_more = true;
    }
    // TODO: no early stopping if the one mutation for the read is
    // to a different value than what read saw in ann_trace
  }
  // DISABLED UNTIL TODO IS SOLVED
  return false;
  if (read_with_more) {
    return false;
  }

  // DEADLOCK-DETECTION-COMPLETENESS
  // If there is a lock that we are about to mutate, no early stopping
  if (!lock_mutations.empty()) {
    return false;
  }
  // Below is hence unreachable
  /*
  std::set<const ZEvent *> maybe_locks;
  for (const auto& l_mut : lock_mutations) {
    assert(is_lock(l_mut.first));
    if (!l_mut.second) {
      // Lock unable to be mutated at this point
      maybe_locks.emplace(l_mut.first);
    } else {
      // Lock with an unlock to mutate to, no early stopping
      return false;
    }
  }
  */

  // If there is a spawn in full and not in part, no early stopping
  for (const auto& cpid_lastspawn : ann_trace.graph().cache().last_spawn) {
    assert(is_spawn(cpid_lastspawn.second));
    assert(ann_trace.graph().has_event(cpid_lastspawn.second));
    assert(ann_trace.po_full().spans_event(cpid_lastspawn.second));
    if (!ann_trace.po_part().spans_event(cpid_lastspawn.second)) {
      return false;
    }
  }

  // DEADLOCK-DETECTION-COMPLETENESS
  // If there is a lock in full and not in part, no early stopping
  for (const auto& cpid_lastlock : ann_trace.graph().cache().last_lock) {
    assert(is_lock(cpid_lastlock.second));
    assert(ann_trace.graph().has_event(cpid_lastlock.second));
    assert(ann_trace.po_full().spans_event(cpid_lastlock.second));
    if (!ann_trace.po_part().spans_event(cpid_lastlock.second)) {
      return false;
    }
  }

  // Consider early stopping

  // Collect threads with a read as the next unannotated event in po_part,
  // such that the read has no negative-allowed mutation even in po_full
  std::set<CPid> blocked_threads;
  for (const ZEvent * read : maybe_reads) {
    assert(is_read(read));
    assert(read_mutations.at(read).empty());
    assert(ann_trace.po_part().spans_event(read));
    assert(!blocked_threads.count(read->cpid()));
    std::set<ZAnn> on_full = ann_trace.graph().mutation_candidates_grouped
    (ann_trace.po_full(), read, ann_trace.negative(), &ann_trace.po_part());
    if (on_full.empty()) {
      // This read has nothing allowed to see in entire full trace
      blocked_threads.insert(read->cpid());
    } else if (on_full.size() > 1) {
      // This read will be able to see different values, no early stopping
      ++early_failed;
      return false;
    }
  }

  // DEADLOCK-DETECTION-COMPLETENESS: this becomes unreachable
  /*
  for (const ZEvent * lock : maybe_locks) {
    assert(is_lock(lock));
    assert(!blocked_threads.count(lock->cpid()));
    const auto& last_unlock = ann_trace.graph().cache().last_unlock;
    if (!last_unlock.count(lock->ml())) {
      blocked_threads.insert(lock->cpid());
      continue;
    }
    assert(last_unlock.count(lock->ml()));
    for (const auto& cpid_lastunlock : last_unlock.at(lock->ml())) {
      assert(is_unlock(cpid_lastunlock.second));
      assert(same_ml(lock, cpid_lastunlock.second));
      assert(ann_trace.po_full().spans_event(cpid_lastunlock.second));
      if (cpid_lastunlock.first == lock->cpid())
        continue;
      if (!ann_trace.po_part().spans_event(cpid_lastunlock.second)) {
        // There will be an unlock to mutate to later on
        // DEADLOCK-DETECTION-COMPLETENESS -- no early stopping
        ++early_failed;
        return false;
      }
    }
    // There will be no unlock to mutate to later on
    blocked_threads.insert(lock->cpid());
  }
  */

  // Blocked threads collected, check all reads of the full trace
  for (const auto& read_localwrite : ann_trace.graph().cache().local_write) {
    const ZEvent * read = read_localwrite.first;
    assert(is_read(read));
    assert(ann_trace.po_full().spans_event(read));
    if (ann_trace.annotation().defines(read)) {
      // This read is already annotated
      assert(ann_trace.po_part().spans_event(read));
      continue;
    }
    if (maybe_reads.count(read)) {
      // This read has already been checked wrt full trace
      continue;
    }
    if (blocked_threads.count(read->cpid())) {
      // This read is preceded by a 'blocked' read/lock, unless that one
      // gets a chance to get 'unblocked', this read will not be reached
      continue;
    }
    assert(read_mutations.count(read) ||
           !ann_trace.po_part().spans_event(read));
    assert(read_mutations.count(read) ||
           !ann_trace.negative().forbids_initial(read));
    std::set<ZAnn> on_full = ann_trace.graph().mutation_candidates_grouped
    (ann_trace.po_full(), read, ann_trace.negative(), &ann_trace.po_part());
    // Separate cases by read in po_part or not
    assert(read_mutations.count(read) || !on_full.empty());
    if (on_full.size() > 1) {
      // This read will be able to see different values, no early stopping
      ++early_failed;
      return false;
    }
    if (!read_mutations.count(read)) {
      // Read outside of po_part with 1 value to see, keep checking
      assert(!ann_trace.po_part().spans_event(read));
      assert(on_full.size() == 1);
      continue;
    }
    // Read on the edge of po_part with one value to see there
    assert(ann_trace.po_part().spans_event(read));
    assert(read_mutations.at(read).size() == 1);
    if (on_full.empty()) {
      // Read has nothing to see in not-in-po_part, keep checking
      continue;
    }
    // One value to see in the not-in-po_part events
    // If that value is different from the one in po_part, no early stopping
    assert(on_full.size() == 1);
    int value_on_part = read_mutations.at(read).begin()->value;
    int value_on_full = on_full.begin()->value;
    if (value_on_part != value_on_full) {
      // This read will be able to see different values, no early stopping
      ++early_failed;
      return false;
    }
  }

  // Early stopping succeeded
  if (blocked_threads.empty()) {
    // We count a full trace since nothing is blocked, i.e.,
    // everything will be negative-allowed to see exactly one value
    ++executed_traces_full;
    // TODO: count full trace only when each of the mutations of
    // read_mutations is realizable
  }
  // TODO: collect backtrack points for all events not in po_part

  ++early_succeeded;
  return true;
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
