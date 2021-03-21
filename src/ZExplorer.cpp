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
static const bool INFO = false;
#include "ZDebug.h"


ZExplorer::ZExplorer(ZBuilderSC& tb)
  : original_tb(&tb), model(MemoryModel::SC) {}


ZExplorer::ZExplorer(ZBuilderTSO& tb)
  : original_tb(&tb), model(MemoryModel::TSO) {}


ZExplorer::ZExplorer(ZBuilderPSO& tb)
  : original_tb(&tb), model(MemoryModel::PSO) {}


void ZExplorer::print_stats() const
{
  std::cout << std::setprecision(2) << std::fixed;
  std::cout << "\n";
  std::cout << "reads_lower_bound_to_check:            " << lin_read_lower_bound << "\n";
  std::cout << "number_of_failed_linearizations:       " << linearization_failed << "\n";
  std::cout << "number_of_cases_below_bound:           " << lin_below_bound << "\n";
  std::cout << "number_of_cases_to_check:              " << lin_goal << "\n";
  std::cout << "number_of_cases_checked:               " << lin_performed << "\n";
  std::cout << "number_of_cases_skipped:               " << lin_skipped << "\n";
  std::cout << "some_algo_spurious_lin:               ";
  for (const int& i : problematic_lin) { std::cout << " " << i; }
  std::cout << "\n";
  assert(lin_performed + lin_below_bound + lin_skipped + linearization_failed == total_lin);
  std::cout << "number_of_events:                     ";
  for (const int& i : no_allevents) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "number_of_reads:                      ";
  for (const int& i : no_reads) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "number_of_writes:                     ";
  for (const int& i : no_writes) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "number_of_threads:                    ";
  for (const int& i : no_threads) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "number_of_variables:                  ";
  for (const int& i : no_variables) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "number_of_closure_rule23_edges:       ";
  for (const int& i : no_closure_rule23_edges) { std::cout << " " << i; }
  std::cout << "\n";
  // Times
  std::cout << "times_rule1_edges:                    ";
  for (const double& d : t_rule1) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "times_closure:                        ";
  for (const double& d : t_closure) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "times_our_yesclosure_yesauxtrace:     ";
  for (const double& d : t_our_yescl_yesaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "times_our_yesclosure_noauxtrace:      ";
  for (const double& d : t_our_yescl_noaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "times_our_noclosure_yesauxtrace:      ";
  for (const double& d : t_our_nocl_yesaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "times_our_noclosure_noauxtrace:       ";
  for (const double& d : t_our_nocl_noaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "times_base_yesclosure_yesauxtrace:    ";
  for (const double& d : t_base_yescl_yesaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "times_base_yesclosure_noauxtrace:     ";
  for (const double& d : t_base_yescl_noaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "times_base_noclosure_yesauxtrace:     ";
  for (const double& d : t_base_nocl_yesaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "times_base_noclosure_noauxtrace:      ";
  for (const double& d : t_base_nocl_noaux) { std::cout << " " << d; }
  std::cout << "\n";
  // branching factors
  /*
  std::cout << "branch_our_yesclosure_yesauxtrace:    ";
  for (const double& d : br_our_yescl_yesaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "branch_our_yesclosure_noauxtrace:     ";
  for (const double& d : br_our_yescl_noaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "branch_our_noclosure_yesauxtrace:     ";
  for (const double& d : br_our_nocl_yesaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "branch_our_noclosure_noauxtrace:      ";
  for (const double& d : br_our_nocl_noaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "branch_base_yesclosure_yesauxtrace:   ";
  for (const double& d : br_base_yescl_yesaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "branch_base_yesclosure_noauxtrace:    ";
  for (const double& d : br_base_yescl_noaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "branch_base_noclosure_yesauxtrace:    ";
  for (const double& d : br_base_nocl_yesaux) { std::cout << " " << d; }
  std::cout << "\n";
  std::cout << "branch_base_noclosure_noauxtrace:     ";
  for (const double& d : br_base_nocl_noaux) { std::cout << " " << d; }
  std::cout << "\n";
  */
  // branching parents and children
  std::cout << "parents_our_yesclosure_yesauxtrace:   ";
  for (const int& i : par_our_yescl_yesaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "children_our_yesclosure_yesauxtrace:  ";
  for (const int& i : ch_our_yescl_yesaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "parents_our_yesclosure_noauxtrace:    ";
  for (const int& i : par_our_yescl_noaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "children_our_yesclosure_noauxtrace:   ";
  for (const int& i : ch_our_yescl_noaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "parents_our_noclosure_yesauxtrace:    ";
  for (const int& i : par_our_nocl_yesaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "children_our_noclosure_yesauxtrace:   ";
  for (const int& i : ch_our_nocl_yesaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "parents_our_noclosure_noauxtrace:     ";
  for (const int& i : par_our_nocl_noaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "children_our_noclosure_noauxtrace:    ";
  for (const int& i : ch_our_nocl_noaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "parents_base_yesclosure_yesauxtrace:  ";
  for (const int& i : par_base_yescl_yesaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "children_base_yesclosure_yesauxtrace: ";
  for (const int& i : ch_base_yescl_yesaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "parents_base_yesclosure_noauxtrace:   ";
  for (const int& i : par_base_yescl_noaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "children_base_yesclosure_noauxtrace:  ";
  for (const int& i : ch_base_yescl_noaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "parents_base_noclosure_yesauxtrace:   ";
  for (const int& i : par_base_nocl_yesaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "children_base_noclosure_yesauxtrace:  ";
  for (const int& i : ch_base_nocl_yesaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "parents_base_noclosure_noauxtrace:    ";
  for (const int& i : par_base_nocl_noaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "children_base_noclosure_noauxtrace:   ";
  for (const int& i : ch_base_nocl_noaux) { std::cout << " " << i; }
  std::cout << "\n";
  std::cout << "\n";
  return;

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


void ZExplorer::dump_schedules() const
{
  llvm::errs() << "SCHEDULES:";
  for (const auto& tauidx_sch : schedules) {
    llvm::errs() << " " << tauidx_sch.first;
  }
  llvm::errs() << "\n";
}


void ZExplorer::maintain_buffers
(BuffersT& buffers, ZEvent * const ev, bool set_up_pointers) const
{
  // Maintain buffers
  if (!buffers.count(ev->cpid().get_proc_seq()))
    buffers.emplace(ev->cpid().get_proc_seq(),
                    std::unordered_map<SymAddrSize, std::list<ZEvent *>>());
  if (isWriteB(ev)) {
    if (!buffers[ev->cpid().get_proc_seq()].count(ev->ml))
      buffers[ev->cpid().get_proc_seq()].emplace(ev->ml, std::list<ZEvent *>());
    buffers[ev->cpid().get_proc_seq()][ev->ml].push_back(ev);
  }
  if (isWriteM(ev)) {
    if (set_up_pointers) {
      // Set up proper WriteB <-> WriteM pointers
      assert(buffers.count(ev->cpid().get_proc_seq()) &&
             buffers.at(ev->cpid().get_proc_seq()).count(ev->ml));
      assert(!buffers.at(ev->cpid().get_proc_seq()).at(ev->ml).empty());
      ZEvent * evB = buffers[ev->cpid().get_proc_seq()][ev->ml].front();
      assert(isWriteB(evB) && sameMl(ev, evB));
      ev->write_other_ptr = evB;
      ev->write_other_trace_id = evB->trace_id();
      evB->write_other_ptr = ev;
      evB->write_other_trace_id = ev->trace_id();
      buffers[ev->cpid().get_proc_seq()][ev->ml].pop_front();
    } else {
      // Assert proper WriteB <-> WriteM pointers
      assert(buffers.count(ev->cpid().get_proc_seq()) &&
             buffers.at(ev->cpid().get_proc_seq()).count(ev->ml));
      assert(!buffers.at(ev->cpid().get_proc_seq()).at(ev->ml).empty());
      ZEvent * evB = buffers[ev->cpid().get_proc_seq()][ev->ml].front();
      assert(isWriteB(evB) && sameMl(ev, evB));
      assert(ev->write_other_ptr == evB);
      assert(ev->write_other_trace_id == evB->trace_id());
      assert(evB->write_other_ptr == ev);
      assert(evB->write_other_trace_id == ev->trace_id());
      buffers[ev->cpid().get_proc_seq()][ev->ml].pop_front();
    }
  }
}


/* *************************** */
/* EXTEND AND EXPLORE          */
/* *************************** */

bool ZExplorer::extend_and_explore
(ZTrace& ann_trace, ZTraceExtension&& ext)
{
  if (lin_performed >= lin_goal) {
    return false;
  }
  start_err("extend...");
  executed_traces++;
  executed_traces_full++;
  assume_blocked_thread += ext.has_assume_blocked_thread;
  executed_traces_full_deadlock += ext.has_deadlocked_thread;
  // Extend ann_trace.exec with ext
  ann_trace.ext_from_id = ext.ext_from_id;
  assert(ann_trace.exec.empty()); // Was removed during get_extension
  ann_trace.exec = std::move(ext.extension);
  assert(ext.extension.empty());
  assert(ann_trace.ext_from_id <= (int) ann_trace.exec.size());
  assert(ann_trace.ext_reads_locks.empty());
  BuffersT buffers; // Thread -> ML -> Pointers to WriteB's
  std::unordered_map<SymAddrSize, int> last_lock; // ML -> ID of last lock to that ML
  readlock_ids.clear();
  previous_lock_id.clear();
  for (int i = 0; i < ann_trace.tau.size(); ++i) {
    ZEvent * ev = ann_trace.tau.at(i).get();
    assert(ev->trace_id() == i);
    // Maintain buffers, Assert proper WriteB <-> WriteM pointers
    maintain_buffers(buffers, ev, false);
    if (isLock(ev)) {
      // Add previous_lock_id
      if (!last_lock.count(ev->ml))
        previous_lock_id.emplace(ev->trace_id(), -1);
      else
        previous_lock_id.emplace(ev->trace_id(), last_lock.at(ev->ml));
      // Maintain last_lock
      assert(!ev->failed_lock);
      if (last_lock.count(ev->ml))
        last_lock.erase(ev->ml);
      last_lock.emplace(ev->ml, ev->trace_id());
    }
    // Assert annotation
    assert(!isRead(ev) || ann_trace.annotation.defines(ev));
    assert(!isLock(ev) || ann_trace.annotation.lock_defines(ev));
    // All readlocks should have readlock_ids
    // Noncommitted readlocks should have schedules
    // Committed readlocks might still have schedules!!!
    if ((isRead(ev) || isLock(ev))) {
      assert(!readlock_ids.count(i));
      readlock_ids.emplace(i);
      if (!ann_trace.committed.count(ev->id())) {
        assert(schedules.count(i));
        assert(failed_schedules.count(i));
        assert(done_schedules.count(i));
      }
    }
  }
  assert(ann_trace.ext_from_id == ann_trace.tau.size());
  // Extend ann_trace.tau/annotation/ext_reads_locks
  for (int i = ann_trace.ext_from_id; i < (int) ann_trace.exec.size(); ++i) {
    ann_trace.tau.push_back(std::unique_ptr<ZEvent>(new ZEvent(
      ann_trace.exec[i], i, true)));
    ZEvent * ev = ann_trace.tau.back().get();
    // Fix the thread_id
    if (ann_trace.proc_seq_to_thread_id.count(ev->cpid().get_proc_seq()))
      ev->_thread_id = ann_trace.proc_seq_to_thread_id[ev->cpid().get_proc_seq()];
    else {
      ev->_thread_id = ann_trace.proc_seq_to_thread_id.size();
      ann_trace.proc_seq_to_thread_id.emplace(
        ev->cpid().get_proc_seq(), ev->thread_id());
    }
    // Maintain buffers, Set up proper WriteB <-> WriteM pointers
    maintain_buffers(buffers, ev, true);
    if (isLock(ev)) {
      // Add previous_lock_id
      if (!last_lock.count(ev->ml))
        previous_lock_id.emplace(ev->trace_id(), -1);
      else
        previous_lock_id.emplace(ev->trace_id(), last_lock.at(ev->ml));
      // Maintain last_lock
      if (!ev->failed_lock) {
        if (last_lock.count(ev->ml))
          last_lock.erase(ev->ml);
        last_lock.emplace(ev->ml, ev->trace_id());
      }
      // Record original_lock
      assert(!original_lock.count(ev->trace_id()));
      original_lock.emplace(ev->trace_id(), ev->id());
    }
    if (isRead(ev) || isLock(ev)) {
      ann_trace.ext_reads_locks.push_back(i);
      assert(!readlock_ids.count(i));
      readlock_ids.emplace(i);
      // Initialize schedules
      assert(!ann_trace.committed.count(ev->id()));
      assert(!schedules.count(i));
      schedules.emplace(i, std::map<ZAnnotation, ZTrace>());
      assert(!failed_schedules.count(i));
      failed_schedules.emplace(i, std::set<ZAnnotation>());
      assert(!done_schedules.count(i));
      done_schedules.emplace(i, std::set<ZAnnotation>());
      // Add observation to annotation
      assert(!ann_trace.annotation.defines(ev) &&
             !ann_trace.annotation.lock_defines(ev));
      if (isRead(ev)) {
        if (ann_trace.exec[i].observed_trace_id == -1)
          ann_trace.annotation.add(ev->id(), ZEventID(true));
        else {
          assert(ann_trace.exec[i].observed_trace_id >= 0 &&
                 ann_trace.exec[i].observed_trace_id < i);
          const ZEvent& obs_ev = ann_trace.exec[ ann_trace.exec[i].observed_trace_id ];
          assert((isWriteB(obs_ev) || isWriteM(obs_ev)) && sameMl(obs_ev, *ev));
          if (isWriteB(obs_ev))
            ann_trace.annotation.add(ev->id(), obs_ev.id());
          else {
            assert(obs_ev.write_other_trace_id < obs_ev.trace_id());
            assert(isWriteB(ann_trace.exec[ obs_ev.write_other_trace_id ]) &&
                   sameMl(*ev, ann_trace.exec[ obs_ev.write_other_trace_id ]));
            ann_trace.annotation.add(ev->id(), ann_trace.exec[ obs_ev.write_other_trace_id ].id());
          }
        }
      } else {
        assert(isLock(ev));
        if (ann_trace.exec[i].observed_trace_id == -1)
          ann_trace.annotation.lock_add(ev->id(), ZEventID(true));
        else {
          assert(ann_trace.exec[i].observed_trace_id >= 0 &&
                 ann_trace.exec[i].observed_trace_id < i);
          const ZEvent& obs_ev = ann_trace.exec[ ann_trace.exec[i].observed_trace_id ];
          assert(isUnlock(obs_ev) && sameMl(obs_ev, *ev));
          ann_trace.annotation.lock_add(ev->id(), obs_ev.id());
        }
      }
    }
  }
  // Explore mutations with extended ann_trace
  if (INFO) { ann_trace.dump(); dump_schedules(); }
  end_err("extend-done");
  return explore(ann_trace);
}


/* *************************** */
/* EXPLORE                     */
/* *************************** */

bool ZExplorer::explore(const ZTrace& ann_trace)
{
  start_err("explore...");
  ZGraph graph(model);
  int mut_start_of_this_explore = mutations_considered;
  for (std::map<int, std::map<ZAnnotation, ZTrace>>::reverse_iterator it = schedules.rbegin();
        it != schedules.rend(); ++it ) {
  //for (const auto& tauidx_sch : schedules) {
    int tauidx = it->first;
    assert(tauidx >= schedules.begin()->first);
    assert(tauidx >= 0 && tauidx < ann_trace.tau.size());
    const ZEvent * ev = ann_trace.tau[tauidx].get();
    assert(ev->trace_id() == tauidx);
    assert(isRead(ev) || isLock(ev));
    if (ann_trace.committed.count(ev->id()))
      continue;
    if (isLock(ev) && previous_lock_id.at(ev->trace_id()) < 0)
      continue;
    if (isLock(ev)) {
      assert(previous_lock_id.at(ev->trace_id()) >= 0);
      if (!original_lock.count(previous_lock_id.at(ev->trace_id()))) {
        assert(ann_trace.committed.count(
               ann_trace.tau.at(previous_lock_id.at(ev->trace_id())).get()->id()));
        continue;
      }
      assert(original_lock.count(previous_lock_id.at(ev->trace_id())));
      if (ev->id() == original_lock.at(previous_lock_id.at(ev->trace_id())))
        continue;
    }
    err_msg(std::string("explore-") + ev->to_string() + "...");
    if (!graph.constructed) {
      clock_t init = std::clock();
      // Construct thread order
      graph.construct(
        ann_trace.tau, ann_trace.tau.size(),
        std::set<int>(), -1);
      // Add reads-from edges
      graph.add_reads_from_edges(ann_trace.annotation);
      time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;
      //if (INFO) graph.dump();
    }
    assert(graph.constructed);
    // Find mutation candidates for read/lock
    bool is_from_extension = ev->trace_id() >= ann_trace.ext_from_id;
    int mutations_only_from_idx = is_from_extension ?
      -1 : ann_trace.ext_from_id;
    std::set<ZEventID> mutations = isRead(ev) ?
      graph.get_read_mutations(ev, ann_trace.annotation.obs(ev), mutations_only_from_idx) :
      graph.get_lock_mutation(ev, previous_lock_id, ann_trace.tau);
    // Perform the mutations
    for (const ZEventID& mutation : mutations) {
      if (lin_performed < lin_goal) {
        mutate(ann_trace, graph, ev, mutation);
      }
    }
  }
  assert(mut_start_of_this_explore <= mutations_considered);
  no_mut_choices += (mut_start_of_this_explore == mutations_considered);
  end_err("explore-done");
  return recur(ann_trace);
}


/* *************************** */
/* MUTATE                      */
/* *************************** */

void ZExplorer::mutate
(const ZTrace& ann_trace, const ZGraph& graph,
 const ZEvent * const readlock, const ZEventID& mutation)
{
  start_err(std::string("mutate-") + readlock->to_string() +
            "   -to-see-   " + mutation.to_string() + "...");
  assert(isRead(readlock) || isLock(readlock));
  assert(!isLock(readlock) ||
         (previous_lock_id.count(readlock->trace_id()) &&
          previous_lock_id.at(readlock->trace_id()) < readlock->trace_id() &&
          previous_lock_id.at(readlock->trace_id()) >= 0));
  int pre_tau_limit = isRead(readlock) ? readlock->trace_id()
                      : previous_lock_id.at(readlock->trace_id());
  assert(pre_tau_limit >= 0);
  auto causes_after = graph.get_causes_after(pre_tau_limit,
    (isRead(readlock) ? mutation : readlock->id()), ann_trace.tau);
  std::set<int>& causes_all_idx = causes_after.first;
  std::set<const ZEvent *>& causes_readslocks = causes_after.second;
  // Key - annotation only for readlock and causes_readslocks
  ZAnnotation mutated_key;
  if (isRead(readlock))
    mutated_key.add(readlock->id(), mutation);
  else {
    assert(isLock(readlock));
    // mutation points to the previous lock, but we need to add to
    // mutated_key the unlock that this previous lock observed
    assert(ann_trace.annotation.lock_defines(mutation));
    const ZEventID& unlock = ann_trace.annotation.lock_obs(mutation);
    mutated_key.lock_add(readlock->id(), unlock);
  }
  for (const auto& cause : causes_readslocks) {
    assert(isRead(cause) || isLock(cause));
    assert(!mutated_key.defines(cause) && !mutated_key.lock_defines(cause));
    if (isRead(cause)) {
      assert(ann_trace.annotation.defines(cause));
      mutated_key.add(cause->id(), ann_trace.annotation.obs(cause));
    } else {
      assert(ann_trace.annotation.lock_defines(cause));
      mutated_key.lock_add(cause->id(), ann_trace.annotation.lock_obs(cause));
    }
  }
  // Stop if such a mutation has already been attempted
  // pre_tau_limit serves as an index to schedules
  assert((isRead(readlock) && pre_tau_limit == readlock->trace_id()) ||
         (isLock(readlock) && pre_tau_limit == previous_lock_id.at(readlock->trace_id())));
  assert(schedules.count(pre_tau_limit) &&
         failed_schedules.count(pre_tau_limit) &&
         done_schedules.count(pre_tau_limit));
  if (failed_schedules[pre_tau_limit].count(mutated_key)) {
    end_err("mutate-done-schedule-already-failed");
    return;
  }
  if (done_schedules[pre_tau_limit].count(mutated_key)) {
    end_err("mutate-done-schedule-already-done");
    return;
  }
  if (schedules[pre_tau_limit].count(mutated_key)) {
    end_err("mutate-done-schedule-already-succeeded");
    return;
  }
  // We attempt the mutation now
  ++mutations_considered;
  clock_t init = std::clock();
  // Construct full annotation
  err_msg("attempt-annotation");
  ZAnnotation mutated_annotation(mutated_key);
  assert((isRead(readlock) && mutated_annotation.defines(readlock)) ||
         (isLock(readlock) && mutated_annotation.lock_defines(readlock)));
  assert((isRead(readlock) && pre_tau_limit == readlock->trace_id()) ||
         (isLock(readlock) && pre_tau_limit == previous_lock_id.at(readlock->trace_id())));
  for (const int tauidx : readlock_ids) {
    assert(tauidx <= pre_tau_limit);
    if (tauidx >= pre_tau_limit) {
      assert(tauidx == pre_tau_limit);
      break;
    }
    assert(tauidx >= 0 && tauidx < (int) ann_trace.tau.size());
    const ZEvent * pre_readlock = ann_trace.tau[tauidx].get();
    assert(isRead(pre_readlock) || isLock(pre_readlock));
    assert(!mutated_annotation.defines(pre_readlock) &&
           !mutated_annotation.lock_defines(pre_readlock));
    if (isRead(pre_readlock)) {
      assert(ann_trace.annotation.defines(pre_readlock->id()));
      mutated_annotation.add(
        pre_readlock->id(), ann_trace.annotation.obs(pre_readlock->id()));
    } else {
      assert(isLock(pre_readlock));
      assert(ann_trace.annotation.lock_defines(pre_readlock->id()));
      mutated_annotation.lock_add(
        pre_readlock->id(), ann_trace.annotation.lock_obs(pre_readlock->id()));
    }
  }
  // Construct the mutation graph
  err_msg("attempt-graph");
  ZGraph mutated_graph(model);
  auto missing_memory_writes = mutated_graph.construct(ann_trace.tau,
    pre_tau_limit, causes_all_idx, readlock->trace_id());
  ZPartialOrder thread_order(mutated_graph.getPo());
  // We add reads-from edges to the thread order
  clock_t r1init = std::clock();
  thread_order.add_reads_from_edges(mutated_annotation);
  double r1time = (double)(clock() - r1init)/CLOCKS_PER_SEC;
  for (int idx : missing_memory_writes) {
    assert(!causes_all_idx.count(idx));
    causes_all_idx.emplace(idx);
  }
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;
  // Close
  err_msg("attempt-closure");
  init = std::clock();
  ZClosure closure(mutated_annotation, mutated_graph, mutated_graph.po);
  bool closed = closure.close();
  double time = (double)(clock() - init)/CLOCKS_PER_SEC;
  time_closure += time;
  if (!closed) {
    ++closure_failed;
    assert(!failed_schedules.at(pre_tau_limit).count(mutated_key));
    failed_schedules[pre_tau_limit].emplace(mutated_key);
    end_err("mutate-done-closure-failed");
    return;
  }
  ++closure_succeeded;
  if (closure.added_edges == 0) {
    closure_no_edge++;
    time_closure_no_edge += time;
  }
  int cur_edges = closure.added_edges;
  closure_edges += closure.added_edges;
  closure_iter += closure.iterations;
  // Linearize
  // if (INFO) { mutated_graph.dump(); mutated_annotation.dump(); }
  err_msg("attempt-linearization");
  ZLinearization linearizer(
    mutated_annotation, mutated_graph.getPo(), ann_trace.exec);
  auto linear = (model != MemoryModel::PSO) ? linearizer.linearizeTSO()
                                            : linearizer.linearizePSO();
  time_linearization += linearizer.elapsed_time;
  ++total_lin;
  total_parents += linearizer.num_parents;
  total_children += linearizer.num_children;
  double cur_br = (linearizer.num_parents==0) ? 1.0 :
    ((double)linearizer.num_children/linearizer.num_parents);
  avg2_branch += ((double)(cur_br - avg2_branch)/total_lin);
  if(cur_br>max_branch) {
    max_branch = cur_br;
  }
  if (linear.empty()) {
    ++linearization_failed;
    assert(!failed_schedules.at(pre_tau_limit).count(mutated_key));
    failed_schedules[pre_tau_limit].emplace(mutated_key);
    end_err("mutate-done-linearization-failed");
    return;
  }
  ++linearization_succeeded;
  //
  //
  //
  //
  //
  if (mutated_annotation.read_size() >= lin_read_lower_bound &&
      (lin_below_bound + lin_skipped) >= (lin_performed * lin_perform_one_per)) {
    // PERFORM LINEARIZATION EXPERIMENTS
    t_rule1.push_back(r1time);
    t_closure.push_back(time);
    t_our_yescl_yesaux.push_back(linearizer.elapsed_time);
    br_our_yescl_yesaux.push_back(cur_br);
    par_our_yescl_yesaux.push_back(linearizer.num_parents);
    ch_our_yescl_yesaux.push_back(linearizer.num_children);
    no_closure_rule23_edges.push_back(cur_edges);
    collect_linearization_stats(linear, mutated_graph);
    linearization_experiments(
      ann_trace, mutated_annotation, mutated_graph.getPo(),
      thread_order);
    //
    // mutated_graph.getPo().dump();
    // thread_order.dump();
    // mutated_annotation.dump();
    //
#ifndef NDEBUG
    ZClosure xxclosure(mutated_annotation, mutated_graph, thread_order);
    bool xxclosed = xxclosure.close();
    assert(xxclosed);
    assert(xxclosure.added_edges == cur_edges);
#endif
  } else {
    if (mutated_annotation.read_size() < lin_read_lower_bound)
      ++lin_below_bound;
    else
      ++lin_skipped;
  }
  //
  //
  //
  //
  //
  assert(total_lin==linearization_succeeded+linearization_failed);
  assert(linearization_respects_ann(linear, mutated_annotation, mutated_graph, ann_trace));
  // Construct tau
  init = std::clock();
  err_msg("attempt-tau");
  std::vector<std::unique_ptr<ZEvent>> mutated_tau;
  BuffersT buffers; // Thread -> ML -> Pointers to WriteB's
  // Add events until pre_tau_limit
  assert((isRead(readlock) && pre_tau_limit == readlock->trace_id()) ||
         (isLock(readlock) && pre_tau_limit == previous_lock_id.at(readlock->trace_id())));
  for (int i = 0; i < pre_tau_limit; ++i) {
    mutated_tau.push_back(std::unique_ptr<ZEvent>(new ZEvent(
      *ann_trace.tau[i].get(), i, true)));
    ZEvent * const ev = mutated_tau.back().get();
    if (isLock(ev)) ev->failed_lock = false;
    else assert(!ev->failed_lock);
    // Maintain buffers, Set up proper WriteB <-> WriteM pointers
    maintain_buffers(buffers, ev, true);
  }
  // Now add the readlock itself
  mutated_tau.push_back(std::unique_ptr<ZEvent>(new ZEvent(
    *ann_trace.tau[readlock->trace_id()].get(), mutated_tau.size(), true)));
  ZEvent * const ev = mutated_tau.back().get();
  if (isLock(ev)) ev->failed_lock = false;
  else assert(!ev->failed_lock);
  // Maintain buffers, Set up proper WriteB <-> WriteM pointers
  maintain_buffers(buffers, ev, true);
  // Now add causes after pre_tau_limit
  for (int cause_idx : causes_all_idx) {
    mutated_tau.push_back(std::unique_ptr<ZEvent>(new ZEvent(
      *ann_trace.tau[cause_idx].get(), mutated_tau.size(), true)));
    ZEvent * const ev = mutated_tau.back().get();
    if (isLock(ev)) ev->failed_lock = false;
    else assert(!ev->failed_lock);
    // Maintain buffers, Set up proper WriteB <-> WriteM pointers
    maintain_buffers(buffers, ev, true);
  }
  // Construct committed reads
  err_msg("attempt-committed");
  std::set<ZEventID> mutated_committed;
  for (const auto& ev_id : ann_trace.committed) {
    assert(graph.hasEvent(ev_id));
    if (mutated_graph.hasEvent(ev_id))
      mutated_committed.emplace(ev_id);
  }
  assert(!ann_trace.committed.count(readlock->id()));
  // Careful: readlock itself does not become
  // committed! That would violate completeness
  for (const auto& cause : causes_readslocks) {
    if (!mutated_committed.count(cause->id()))
      mutated_committed.emplace(cause->id());
  }
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;
  // Add successful schedule
  assert((isRead(readlock) && pre_tau_limit == readlock->trace_id()) ||
         (isLock(readlock) && pre_tau_limit == previous_lock_id.at(readlock->trace_id())));
  assert(schedules.count(pre_tau_limit) &&
         failed_schedules.count(pre_tau_limit) &&
         done_schedules.count(pre_tau_limit));
  assert(!schedules.at(pre_tau_limit).count(mutated_key) &&
         !failed_schedules.at(pre_tau_limit).count(mutated_key) &&
         !done_schedules.at(pre_tau_limit).count(mutated_key));
  schedules[pre_tau_limit].emplace(mutated_key,
    ZTrace(&ann_trace, std::move(linear), std::move(mutated_tau),
           std::move(mutated_annotation), std::move(mutated_committed),
           readlock->id()));
  end_err(std::string("mutate-done-added-to-") + std::to_string(pre_tau_limit));
}


/* *************************** */
/* RECUR                       */
/* *************************** */

bool ZExplorer::recur(const ZTrace& ann_trace)
{
  start_err("recur...");
  // Iterate through reads/locks of extension in reverse order
  for (int i = ann_trace.ext_reads_locks.size() - 1; i >= 0; --i) {
    int idx = ann_trace.ext_reads_locks[i];
    assert(idx < ann_trace.tau.size());
    assert(isRead(ann_trace.tau.at(idx).get()) ||
           isLock(ann_trace.tau.at(idx).get()));
    start_err(std::string("recur-on-") + std::to_string(i) +
            "-(idx-is-" + std::to_string(idx) + ")...");
    assert(schedules.count(idx) && failed_schedules.count(idx) &&
           done_schedules.count(idx));
    while (!schedules.at(idx).empty()) {
      // Move the schedules out of the map
      std::map<ZAnnotation, ZTrace> sch = std::move(schedules.at(idx));
      assert(schedules.at(idx).empty());
      for (const auto& key_trace : sch) {
        assert(!done_schedules.at(idx).count(key_trace.first));
        done_schedules.at(idx).emplace(key_trace.first);
      }
      // Recur on each schedule
      for (auto& key_trace : sch) {
        start_err(std::string("schedule-on-") + std::to_string(i) +
                "-(idx-is-" + std::to_string(idx) + ")...");
        if (lin_performed >= lin_goal) {
          return false;
        }
        if (get_extension(key_trace.second)) {
          return true;
        }
        end_err(std::string("schedule-on-") + std::to_string(i) +
                "-(idx-is-" + std::to_string(idx) + ")-done");
      }
      // While exploring above schedules, additional novel schedules
      // for this same readlock may get introduced. Therefore we
      // check the schedules set again for these novel schedules.
    }
    // We erase all schedules only now
    assert(schedules.count(idx) && schedules.at(idx).empty());
    schedules.erase(idx);
    assert(!schedules.count(idx));
    assert(failed_schedules.count(idx));
    failed_schedules.erase(idx);
    assert(!failed_schedules.count(idx));
    assert(done_schedules.count(idx));
    done_schedules.erase(idx);
    assert(!done_schedules.count(idx));
    // Erase original_lock
    if (isLock(ann_trace.tau.at(idx).get())) {
      assert(original_lock.count(idx));
      original_lock.erase(idx);
    } else
      assert(!original_lock.count(idx));
    end_err(std::string("recur-on-") + std::to_string(i) +
            "-(idx-is-" + std::to_string(idx) + ")-done");
  }
  end_err("recur-done");
  return false;
}


/* *************************** */
/* GET EXTENSION               */
/* *************************** */

bool ZExplorer::get_extension(ZTrace& ann_trace)
{
  start_err("get_extension...");
  assert(original_tb);
  interpreter_used++;
  ZTraceExtension ext;
  //if (INFO) dump_trace(ann_trace.exec);
  if (model == MemoryModel::SC) {
    clock_t init = std::clock();
    ZBuilderSC TB(*(original_tb->config), original_tb->M,
                  std::move(ann_trace.exec));
    ext = TB.extendGivenTrace();
    time_interpreter += (double)(clock() - init)/CLOCKS_PER_SEC;
    // Error
    if (TB.has_error()) {
      original_tb->error_trace = TB.get_trace();
      end_err("get_extension-error-SC");
      return true;
    }
  } else if (model == MemoryModel::TSO) {
    clock_t init = std::clock();
    ZBuilderTSO TB(*(original_tb->config), original_tb->M,
                   std::move(ann_trace.exec));
    ext = TB.extendGivenTrace();
    time_interpreter += (double)(clock() - init)/CLOCKS_PER_SEC;
    // Error
    if (TB.has_error()) {
      original_tb->error_trace = TB.get_trace();
      end_err("get_extension-error-TSO");
      return true;
    }
  } else {
    assert(model == MemoryModel::PSO);
    clock_t init = std::clock();
    ZBuilderPSO TB(*(original_tb->config), original_tb->M,
                   std::move(ann_trace.exec));
    ext = TB.extendGivenTrace();
    time_interpreter += (double)(clock() - init)/CLOCKS_PER_SEC;
    // Error
    if (TB.has_error()) {
      original_tb->error_trace = TB.get_trace();
      end_err("get_extension-error-PSO");
      return true;
    }
  }
  // Extend with the extension and explore
  end_err("get_extension-done");
  //if (INFO) ext.dump();
  return extend_and_explore(ann_trace, std::move(ext));
}


/* *************************** */
/* LINEARIZATION EXPERIMENTS   */
/* *************************** */

void ZExplorer::collect_linearization_stats(
const std::vector<ZEvent>& linearization, const ZGraph& graph)
{
  ++lin_performed;
  if (lin_performed == 1) {
    std::cout << "There are large-enough linearization instances to consider\n";
  }
  no_allevents.push_back(linearization.size());
  no_threads.push_back(int(graph.number_of_threads()));
  std::unordered_set<SymAddrSize> variables;
  int reads = 0;
  int writes = 0;
  for (const ZEvent& ev : linearization) {
    if (isRead(ev)) {
      ++reads;
      if (!variables.count(ev.ml)) {
        variables.insert(ev.ml);
      }
    } else if (isWriteM(ev)) {
      ++writes;
      if (!variables.count(ev.ml)) {
        variables.insert(ev.ml);
      }
    }
  }
  no_reads.push_back(reads);
  no_writes.push_back(writes);
  no_variables.push_back(variables.size());
}


void ZExplorer::linearization_experiments(
const ZTrace& ann_trace, const ZAnnotation& annotation,
const ZPartialOrder& closed_po, const ZPartialOrder& thread_order)
{
  for (int i = 0; i < ann_trace.exec.size(); ++i) {
    ann_trace.exec.at(i).exec_trace_id = i;
    int thi = ann_trace.exec.at(i).thread_id();
    int aui = ann_trace.exec.at(i).aux_id();
    int evi = ann_trace.exec.at(i).event_id();
    if (closed_po.graph.hasThreadAux(thi, aui) &&
        closed_po.graph(thi, aui).size() > evi) {
      closed_po.graph.getEvent(thi, aui, evi)->exec_trace_id = i;
    }
  }

  // OUR, YESclosure, NOauxtrace
  {
  std::vector<ZEvent> emptyaux;
  assert(emptyaux.empty());
  ZLinearization linearizer(
    annotation, closed_po, emptyaux);
  std::vector<ZEvent> linear = (model != MemoryModel::PSO) ?
    linearizer.linearizeTSO() : linearizer.linearizePSO();
  if (!linear.empty()) {
    assert(!linearizer.exceeded_limit);
    assert(linearizer.elapsed_time <= linearizer.time_limit);
    t_our_yescl_noaux.push_back(linearizer.elapsed_time);
    bool respects = linearization_respects_ann(linear, annotation, closed_po.graph, ann_trace);
    if (!respects && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(respects);
  } else {
    if (!linearizer.exceeded_limit && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(linearizer.exceeded_limit);
    assert(linearizer.elapsed_time > linearizer.time_limit);
    t_our_yescl_noaux.push_back(linearizer.time_limit);
  }
  double cur_br = (linearizer.num_parents==0) ? 1.0 :
    ((double)linearizer.num_children/linearizer.num_parents);
  br_our_yescl_noaux.push_back(cur_br);
  par_our_yescl_noaux.push_back(linearizer.num_parents);
  ch_our_yescl_noaux.push_back(linearizer.num_children);
  }

  // OUR, NOclosure, YESauxtrace
  {
  ZLinNoclosure linearizer(
    annotation, thread_order, ann_trace.exec);
  std::vector<ZEvent> linear = (model != MemoryModel::PSO) ?
    linearizer.linearizeTSO() : linearizer.linearizePSO();
  if (!linear.empty()) {
    assert(!linearizer.exceeded_limit);
    assert(linearizer.elapsed_time <= linearizer.time_limit);
    t_our_nocl_yesaux.push_back(linearizer.elapsed_time);
    bool respects = linearization_respects_ann(linear, annotation, closed_po.graph, ann_trace);
    if (!respects && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(respects);
  } else {
    if (!linearizer.exceeded_limit && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(linearizer.exceeded_limit);
    assert(linearizer.elapsed_time > linearizer.time_limit);
    t_our_nocl_yesaux.push_back(linearizer.time_limit);
  }
  double cur_br = (linearizer.num_parents==0) ? 1.0 :
    ((double)linearizer.num_children/linearizer.num_parents);
  br_our_nocl_yesaux.push_back(cur_br);
  par_our_nocl_yesaux.push_back(linearizer.num_parents);
  ch_our_nocl_yesaux.push_back(linearizer.num_children);
  }

  // OUR, NOclosure, NOauxtrace
  {
  std::vector<ZEvent> emptyaux;
  assert(emptyaux.empty());
  ZLinNoclosure linearizer(
    annotation, thread_order, emptyaux);
  std::vector<ZEvent> linear = (model != MemoryModel::PSO) ?
    linearizer.linearizeTSO() : linearizer.linearizePSO();
  if (!linear.empty()) {
    assert(!linearizer.exceeded_limit);
    assert(linearizer.elapsed_time <= linearizer.time_limit);
    t_our_nocl_noaux.push_back(linearizer.elapsed_time);
    bool respects = linearization_respects_ann(linear, annotation, closed_po.graph, ann_trace);
    if (!respects && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(respects);
  } else {
    if (!linearizer.exceeded_limit && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(linearizer.exceeded_limit);
    assert(linearizer.elapsed_time > linearizer.time_limit);
    t_our_nocl_noaux.push_back(linearizer.time_limit);
  }
  double cur_br = (linearizer.num_parents==0) ? 1.0 :
    ((double)linearizer.num_children/linearizer.num_parents);
  br_our_nocl_noaux.push_back(cur_br);
  par_our_nocl_noaux.push_back(linearizer.num_parents);
  ch_our_nocl_noaux.push_back(linearizer.num_children);
  }

  //////////////////////////////////////////////////////////////////

  // BASE, YESclosure, YESauxtrace
  {
  ZLinNaive linearizer(
    annotation, closed_po, ann_trace.exec);
  std::vector<ZEvent> linear = linearizer.linearize();
  if (!linear.empty()) {
    assert(!linearizer.exceeded_limit);
    assert(linearizer.elapsed_time <= linearizer.time_limit);
    t_base_yescl_yesaux.push_back(linearizer.elapsed_time);
    bool respects = linearization_respects_ann(linear, annotation, closed_po.graph, ann_trace);
    if (!respects && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(respects);
  } else {
    if (!linearizer.exceeded_limit && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(linearizer.exceeded_limit);
    assert(linearizer.elapsed_time > linearizer.time_limit);
    t_base_yescl_yesaux.push_back(linearizer.time_limit);
  }
  double cur_br = (linearizer.num_parents==0) ? 1.0 :
    ((double)linearizer.num_children/linearizer.num_parents);
  br_base_yescl_yesaux.push_back(cur_br);
  par_base_yescl_yesaux.push_back(linearizer.num_parents);
  ch_base_yescl_yesaux.push_back(linearizer.num_children);
  }

  // BASE, YESclosure, NOauxtrace
  {
  std::vector<ZEvent> emptyaux;
  assert(emptyaux.empty());
  ZLinNaive linearizer(
    annotation, closed_po, emptyaux);
  std::vector<ZEvent> linear = linearizer.linearize();
  if (!linear.empty()) {
    assert(!linearizer.exceeded_limit);
    assert(linearizer.elapsed_time <= linearizer.time_limit);
    t_base_yescl_noaux.push_back(linearizer.elapsed_time);
    bool respects = linearization_respects_ann(linear, annotation, closed_po.graph, ann_trace);
    if (!respects && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(respects);
  } else {
    if (!linearizer.exceeded_limit && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(linearizer.exceeded_limit);
    assert(linearizer.elapsed_time > linearizer.time_limit);
    t_base_yescl_noaux.push_back(linearizer.time_limit);
  }
  double cur_br = (linearizer.num_parents==0) ? 1.0 :
    ((double)linearizer.num_children/linearizer.num_parents);
  br_base_yescl_noaux.push_back(cur_br);
  par_base_yescl_noaux.push_back(linearizer.num_parents);
  ch_base_yescl_noaux.push_back(linearizer.num_children);
  }

  // BASE, NOclosure, YESauxtrace
  {
  ZLinNaive linearizer(
    annotation, thread_order, ann_trace.exec);
  std::vector<ZEvent> linear = linearizer.linearize();
  if (!linear.empty()) {
    assert(!linearizer.exceeded_limit);
    assert(linearizer.elapsed_time <= linearizer.time_limit);
    t_base_nocl_yesaux.push_back(linearizer.elapsed_time);
    bool respects = linearization_respects_ann(linear, annotation, closed_po.graph, ann_trace);
    if (!respects && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(respects);
  } else {
    if (!linearizer.exceeded_limit && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(linearizer.exceeded_limit);
    assert(linearizer.elapsed_time > linearizer.time_limit);
    t_base_nocl_yesaux.push_back(linearizer.time_limit);
  }
  double cur_br = (linearizer.num_parents==0) ? 1.0 :
    ((double)linearizer.num_children/linearizer.num_parents);
  br_base_nocl_yesaux.push_back(cur_br);
  par_base_nocl_yesaux.push_back(linearizer.num_parents);
  ch_base_nocl_yesaux.push_back(linearizer.num_children);
  }

  // BASE, NOclosure, NOauxtrace
  {
  std::vector<ZEvent> emptyaux;
  assert(emptyaux.empty());
  ZLinNaive linearizer(
    annotation, thread_order, emptyaux);
  std::vector<ZEvent> linear = linearizer.linearize();
  if (!linear.empty()) {
    assert(!linearizer.exceeded_limit);
    assert(linearizer.elapsed_time <= linearizer.time_limit);
    t_base_nocl_noaux.push_back(linearizer.elapsed_time);
    bool respects = linearization_respects_ann(linear, annotation, closed_po.graph, ann_trace);
    if (!respects && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(respects);
  } else {
    if (!linearizer.exceeded_limit && problematic_lin.back() != lin_performed - 1) { problematic_lin.push_back(lin_performed - 1); }
    assert(linearizer.exceeded_limit);
    assert(linearizer.elapsed_time > linearizer.time_limit);
    t_base_nocl_noaux.push_back(linearizer.time_limit);
  }
  double cur_br = (linearizer.num_parents==0) ? 1.0 :
    ((double)linearizer.num_children/linearizer.num_parents);
  br_base_nocl_noaux.push_back(cur_br);
  par_base_nocl_noaux.push_back(linearizer.num_parents);
  ch_base_nocl_noaux.push_back(linearizer.num_children);
  }
}


/* *************************** */
/* RESPECTS ANNOTATION         */
/* *************************** */

bool ZExplorer::linearization_respects_ann
(const std::vector<ZEvent>& trace,
 const ZAnnotation& annotation,
 const ZGraph& graph,
 const ZTrace& parent_trace) const
{
  assert(!trace.empty());
  // Trace index of the last write happening in the given location
  std::unordered_map<SymAddrSize, int> lastWrite;
  // Trace index of the last lock + whether given ML is currently locked
  std::unordered_map<SymAddrSize, int> lastUnlock;
  std::unordered_set<SymAddrSize> locked;
  // ThreadID -> ML -> Writes in store queue of thr for ml
  std::unordered_map
    <int, std::unordered_map
     <SymAddrSize, std::list<int>>> storeQueue;
  // ThreadID -> order in which enqueued
  std::unordered_map<int, std::vector<int>> queue_enqorder;
  // ThreadID -> order in which dequeued
  std::unordered_map<int, std::vector<int>> queue_deqorder;
  // Real observation in trace
  std::unordered_map<int, int> realObs;
  // BufferWrite -> pos
  std::unordered_map<ZEventID, unsigned> bw_pos;
  // BufferWrite -> MemoryWrite
  std::unordered_map<int, int> buf_mem;
  // Number of annotated reads/locks
  int annotated_readlocks = 0;

  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent * ev = &(trace.at(i));
    assert(graph.proc_seq_to_thread_id.count(ev->cpid().get_proc_seq()));
    assert(ev->thread_id() == graph.proc_seq_to_thread_id.at(ev->cpid().get_proc_seq()));

    if (isWriteB(ev)) {
      // bw_pos
      bw_pos.emplace(ev->id(), i);
      // storeQueue
      if (!storeQueue.count(ev->thread_id()))
        storeQueue.emplace
          (ev->thread_id(), std::unordered_map<SymAddrSize, std::list<int>>());
      if (!storeQueue.at(ev->thread_id()).count(ev->ml))
        storeQueue.at(ev->thread_id()).emplace
          (ev->ml, std::list<int>());
      storeQueue.at(ev->thread_id()).at(ev->ml).push_back(i);
      // queue_enqorder
      if (!queue_enqorder.count(ev->thread_id())) {
        queue_enqorder.emplace(ev->thread_id(), std::vector<int>());
      }
      queue_enqorder.at(ev->thread_id()).push_back(i);
    }

    if (isWriteM(ev)) {
      // lastWrite
      lastWrite[ev->ml] = i;
      // storeQueue
      assert(storeQueue.count(ev->thread_id()) && "Maybe writeM before writeB?");
      assert(storeQueue.at(ev->thread_id()).count(ev->ml) && "Maybe writeM before writeB?");
      assert(!storeQueue.at(ev->thread_id()).at(ev->ml).empty());
      buf_mem.emplace(storeQueue.at(ev->thread_id()).at(ev->ml).front(), i);
      //ev->write_other_trace_id = storeQueue.at(ev->thread_id()).at(ev->ml).front(); ////////////////////
      storeQueue.at(ev->thread_id()).at(ev->ml).pop_front();
      // queue_deqorder
      if (!queue_deqorder.count(ev->thread_id())) {
        queue_deqorder.emplace(ev->thread_id(), std::vector<int>());
      }
      queue_deqorder.at(ev->thread_id()).push_back(i);
      assert(queue_enqorder.count(ev->thread_id()) && "Maybe writeM before writeB?");
      assert(queue_enqorder.at(ev->thread_id()).size() >= queue_deqorder.at(ev->thread_id()).size());
      int last = queue_deqorder.at(ev->thread_id()).size() - 1;
      assert(last >= 0);
      const ZEvent& bw = trace.at(queue_enqorder.at(ev->thread_id()).at(last));
      assert(isWriteB(bw));
      if (model != MemoryModel::PSO &&
          ((last >= 1 && ev->event_id() <= trace.at(queue_deqorder.at(ev->thread_id()).at(last - 1)).event_id()) ||
           ev->ml != bw.ml)) {
        llvm::errs() << "Partial order\n";
        graph.dump();
        llvm::errs() << "Linearization\n";
        dump_trace(trace);
        llvm::errs() << "Annotation that should be respected\n";
        annotation.dump();
        llvm::errs() << "Memory-write dequeued earlier than another memory-write added to the buffer beforehand\n";
        ev->dump();
        return false;
      }
    }

    if (isRead(ev)) {
      assert(annotation.defines(ev->id()));
      ++annotated_readlocks;
      if (storeQueue.count(ev->thread_id()) &&
          storeQueue.at(ev->thread_id()).count(ev->ml) &&
          !storeQueue.at(ev->thread_id()).at(ev->ml).empty()) {
        realObs.emplace(i, storeQueue.at(ev->thread_id()).at(ev->ml).back());
        //ev->observed_trace_id = storeQueue.at(ev->thread_id()).at(ev->ml).back(); ////////////////////
        //ev->value = trace.at(ev->observed_trace_id).value; ////////////////////
      }
      else if (lastWrite.count(ev->ml)) {
        realObs.emplace(i, lastWrite.at(ev->ml));
        //ev->observed_trace_id = lastWrite.at(ev->ml); ////////////////////
        //ev->value = trace.at(ev->observed_trace_id).value; ////////////////////
      }
      else {
        realObs.emplace(i, -1);
        //ev->observed_trace_id = -1; ////////////////////
        //ev->value = 0; ////////////////////
      }
    }

    if (isLock(ev)) {
      assert(annotation.lock_defines(ev->id()));
      ++annotated_readlocks;
      if (locked.count(ev->ml)) {
        llvm::errs() << "Partial order\n";
        graph.dump();
        llvm::errs() << "Linearization\n";
        dump_trace(trace);
        llvm::errs() << "Annotation that should be respected\n";
        annotation.dump();
        llvm::errs() << "Lock trying to acquire already locked mutex\n";
        ev->dump();
        return false;
      }
      assert(!locked.count(ev->ml));
      locked.insert(ev->ml);
      // Check whether the lock is consistent with the annotation
      assert(annotation.lock_defines(ev));
      const auto& an = annotation.lock_obs(ev);
      if (!lastUnlock.count(ev->ml)) {
        assert(an == ZEventID(true)); // initial
        //ev->observed_trace_id = -1; ////////////////////
      }
      else {
        assert(an == trace.at(lastUnlock.at(ev->ml)).id());
        //ev->observed_trace_id = lastUnlock.at(ev->ml); ////////////////////
      }
    }
    if (isUnlock(ev)) {
      assert(locked.count(ev->ml));
      locked.erase(ev->ml);
      lastUnlock[ev->ml] = i;
    }
  }
  if (annotated_readlocks != annotation.size()) {
    llvm::errs() << "Partial order\n";
    graph.dump();
    llvm::errs() << "Linearization\n";
    dump_trace(trace);
    llvm::errs() << "Annotation that should be respected\n";
    annotation.dump();
    llvm::errs() << "Annotation has " << annotation.size() <<
      "read/locks, linearization has " << annotated_readlocks << "\n";
    return false;
  }
  for (const auto& bwid_bwpos : bw_pos) {
    if (!buf_mem.count(bwid_bwpos.second)) {
      const ZEvent& bw = trace.at(bwid_bwpos.second);
      assert(isWriteB(bw));
      llvm::errs() << "Partial order\n";
      graph.dump();
      llvm::errs() << "Linearization\n";
      dump_trace(trace);
      llvm::errs() << "Annotation that should be respected\n";
      annotation.dump();
      llvm::errs() << "Buffer-write has no memory-write in the trace";
      bw.dump();
      return false;
    }
  }

  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = &(trace.at(i));
    if (isRead(ev)) {
      assert(annotation.defines(ev->id()));
      const ZEventID& obs = annotation.obs(ev->id());
      const ZEvent *obsB = graph.initial();
      if (obs.event_id() >= 0) {
        assert(bw_pos.count(obs));
        obsB = &(trace.at(bw_pos.at(obs)));
        assert(isWriteB(obsB));
      }
      const ZEvent *obsM = (isInitial(obsB))
        ? graph.initial() : &(trace.at(buf_mem.at(bw_pos.at(obs))));;
      assert(obsB->value == obsM->value);
      assert(isInitial(obsB) || (isWriteB(obsB) && isWriteM(obsM) &&
                                 sameMl(obsB, obsM) && sameMl(ev, obsB)));
      const ZEvent *realObservation = (realObs.at(i) == -1)
        ? graph.initial() : &(trace.at(realObs.at(i)));
      if (*realObservation != *obsB && *realObservation != *obsM) {
        llvm::errs() << "Partial order\n";
        graph.dump();
        llvm::errs() << "Linearization\n";
        dump_trace(trace);
        llvm::errs() << "Annotation that should be respected\n";
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

bool ZExplorer::respectsAnnotation
(const std::vector<ZEvent>& trace,
 const ZAnnotation& annotation,
 const ZPartialOrder& mutatedPO,
 const ZTrace& parentTrace) const
{
  assert(!trace.empty());
  std::unordered_map<ZEventID, unsigned> bw_pos;
  for (const ZEvent& evref: trace) {
    const ZEvent *ev = &evref;
    // ev->thread_id() not set, have to get it
    unsigned thrid = parentTrace.graph.getThreadIDnoAdd(ev);
    ev->_thread_id = thrid;
    if (isWriteB(ev) && thrid != 1337)
      bw_pos.emplace(ev->id(), ev->trace_id());
  }
  for (const ZEvent& evref : trace) {
    const ZEvent *ev = &evref;
    if (isRead(ev)) {
      if (annotation.defines(ev->id())) {
        const ZEventID& obs = annotation.obs(ev->id());
        const ZEvent *obsB = parentTrace.graph.initial();
        if (obs.event_id() >= 0) {
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
          dump_trace(trace);
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

*/
