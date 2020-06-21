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

static const bool DEBUG = true;
static const bool INFO = true;
#include "ZDebug.h"


ZExplorer::ZExplorer(ZBuilderSC& tb)
  : original_tb(&tb), model(MemoryModel::SC) {}


ZExplorer::ZExplorer(ZBuilderTSO& tb)
  : original_tb(&tb), model(MemoryModel::TSO) {}


ZExplorer::ZExplorer(ZBuilderPSO& tb)
  : original_tb(&tb), model(MemoryModel::PSO) {}


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
/* EXTEND AND EXPLORE          */
/* *************************** */

bool ZExplorer::extend_and_explore
(ZTrace& ann_trace, ZTraceExtension&& ext)
{
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
  // Thread -> ML -> Pointers to WriteB's
  std::unordered_map<
    std::vector<int>, std::unordered_map<SymAddrSize, std::list<ZEvent *>>> buffers;
  // Extend ann_trace.tau/annotation/ext_reads_locks
  if (INFO) dump_trace(ann_trace.exec);
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
    if (isRead(ev) || isLock(ev)) {
      ann_trace.ext_reads_locks.push_back(i);
      // Initialize schedules
      assert(!schedules.count(i));
      schedules.emplace(i, std::map<ZAnnotation, ZTrace>());
      assert(!failed_schedules.count(i));
      failed_schedules.emplace(i, std::set<ZAnnotation>());
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
    // Maintain buffers
    if (!buffers.count(ev->cpid().get_proc_seq()))
      buffers.emplace(ev->cpid().get_proc_seq(),
                      std::unordered_map<SymAddrSize, std::list<ZEvent *>>());
    if (isWriteB(ev)) {
      if (!buffers[ev->cpid().get_proc_seq()].count(ev->ml))
        buffers[ev->cpid().get_proc_seq()].emplace(ev->ml, std::list<ZEvent *>());
      buffers[ev->cpid().get_proc_seq()][ev->ml].push_back(ev);
    }
    // Set up proper WriteB <-> WriteM pointers
    if (isWriteM(ev)) {
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
    }
  }
  // Explore mutations with extended ann_trace
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
  for (const auto& tauidx_sch : schedules) {
    assert(tauidx_sch.first >= 0 && tauidx_sch.first < ann_trace.tau.size());
    const ZEvent * ev = ann_trace.tau[tauidx_sch.first].get();
    assert(ev->trace_id() == tauidx_sch.first);
    assert(isRead(ev) || isLock(ev));
    if (ann_trace.committed.count(ev->id()))
      continue;
    err_msg(std::string("explore-") + ev->to_string() + "...");
    if (!graph.constructed) {
      clock_t init = std::clock();
      // Construct thread order
      graph.construct(
        ann_trace.tau, ann_trace.tau.size(), std::set<int>());
      // Add reads-from edges
      graph.add_reads_from_edges(ann_trace.annotation);
      time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;
    }
    assert(graph.constructed);
    // Find mutation candidates for read/lock
    const ZEventID& base_obs = isRead(ev) ?
      ann_trace.annotation.obs(ev): ann_trace.annotation.lock_obs(ev);
    bool is_from_extension = ev->trace_id() >= ann_trace.ext_from_id;
    int mutations_only_from_idx = is_from_extension ?
      -1 : ann_trace.ext_from_id;
    std::set<ZEventID> mutations = graph.get_mutations(
      ev, base_obs, mutations_only_from_idx);
    // Perform the mutations
    for (const ZEventID& mutation : mutations)
      mutate(ann_trace, graph, ev, mutation);
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
  start_err(std::string("mutate") + readlock->to_string() +
            "   -to-see-   " + mutation.to_string() + "...");
  assert(isRead(readlock) || isLock(readlock));
  auto causes_after = graph.get_causes_after(
    readlock, mutation, ann_trace.tau);
  std::set<int>& causes_all_idx = causes_after.first;
  std::set<const ZEvent *>& causes_readslocks = causes_after.second;
  // Key - annotation only for readlock and causes_readslocks
  ZAnnotation mutated_key;
  if (isRead(readlock)) mutated_key.add(readlock->id(), mutation);
  else mutated_key.lock_add(readlock->id(), mutation);
  for (const auto& cause : causes_readslocks) {
    assert(isRead(cause) || isLock(cause));
    assert(!mutated_key.defines(cause) && !mutated_key.lock_defines(cause));
    if (isRead(cause)) {
      assert(ann_trace.annotation.defines(cause));
      mutated_key.add(cause->id(), ann_trace.annotation.obs(cause));
    } else {
      assert(ann_trace.annotation.lock_defines(cause));
      mutated_key.lock_add(cause->id(), ann_trace.annotation.obs(cause));
    }
  }
  // Stop if such a mutation has already been attempted
  assert(schedules.count(readlock->trace_id()) &&
         failed_schedules.count(readlock->trace_id()));
  if (failed_schedules[readlock->trace_id()].count(mutated_key)) {
    end_err("mutate-done-not-added-already-failed");
    return;
  }
  if (schedules[readlock->trace_id()].count(mutated_key)) {
    end_err("mutate-done-not-added-already-succeeded");
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
  for (const auto& tauidx_sch : schedules) {
    int tauidx = tauidx_sch.first;
    assert(tauidx <= readlock->trace_id());
    if (tauidx == readlock->trace_id())
      break;
    assert(tauidx >= 0 && tauidx < (int) ann_trace.tau.size());
    const ZEvent * pre_readlock = ann_trace.tau[tauidx].get();
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
  mutated_graph.construct(
    ann_trace.tau, readlock->trace_id(), causes_all_idx);
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;
  // Close
  err_msg("attempt-closure");
  init = std::clock();
  ZClosure closure(mutated_annotation, mutated_graph);
  bool closed = closure.close();
  double time = (double)(clock() - init)/CLOCKS_PER_SEC;
  time_closure += time;
  if (!closed) {
    ++closure_failed;
    assert(!failed_schedules.at(readlock->trace_id()).count(mutated_key));
    failed_schedules[readlock->trace_id()].emplace(mutated_key);
    end_err("mutate-done-closure-failed");
    return;
  }
  ++closure_succeeded;
  if (closure.added_edges == 0) {
    closure_no_edge++;
    time_closure_no_edge += time;
  }
  closure_edges += closure.added_edges;
  closure_iter += closure.iterations;
  // Linearize
  if (INFO) { mutated_graph.dump(); mutated_annotation.dump(); }
  err_msg("attempt-linearization");
  init = std::clock();
  ZLinearization linearizer(
    mutated_annotation, mutated_graph.getPo(), ann_trace.exec);
  auto linear = (model != MemoryModel::PSO) ? linearizer.linearizeTSO()
                                            : linearizer.linearizePSO();
  time_linearization += (double)(clock() - init)/CLOCKS_PER_SEC;
  total_parents += linearizer.num_parents;
  total_children += linearizer.num_children;
  if (linear.empty()) {
    assert(!failed_schedules.at(readlock->trace_id()).count(mutated_key));
    failed_schedules[readlock->trace_id()].emplace(mutated_key);
    end_err("mutate-done-linearization-failed");
    return;
  }
  // TODO  assert(linearization_respects_ann(linear, mutated_annotation, mutated_graph, ann_trace));  TODO
  // Construct tau
  init = std::clock();
  err_msg("attempt-tau");
  std::vector<std::unique_ptr<ZEvent>> mutated_tau;
  // Thread -> ML -> Pointers to WriteB's
  std::unordered_map<
    std::vector<int>, std::unordered_map<SymAddrSize, std::list<ZEvent *>>> buffers;
  for (int i = 0; i <= readlock->trace_id(); ++i) {
    mutated_tau.push_back(std::unique_ptr<ZEvent>(new ZEvent(
      *ann_trace.tau[i].get(), i, true)));
    ZEvent * const ev = mutated_tau.back().get();
    // Maintain buffers
    if (!buffers.count(ev->cpid().get_proc_seq()))
      buffers.emplace(ev->cpid().get_proc_seq(),
                      std::unordered_map<SymAddrSize, std::list<ZEvent *>>());
    if (isWriteB(ev)) {
      if (!buffers[ev->cpid().get_proc_seq()].count(ev->ml))
        buffers[ev->cpid().get_proc_seq()].emplace(ev->ml, std::list<ZEvent *>());
      buffers[ev->cpid().get_proc_seq()][ev->ml].push_back(ev);
    }
    // Set up proper WriteB <-> WriteM pointers
    if (isWriteM(ev)) {
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
    }
  }
  for (int cause_idx : causes_all_idx) {
    mutated_tau.push_back(std::unique_ptr<ZEvent>(new ZEvent(
      *ann_trace.tau[cause_idx].get(), mutated_tau.size(), true)));
    ZEvent * const ev = mutated_tau.back().get();
    // Maintain buffers
    if (!buffers.count(ev->cpid().get_proc_seq()))
      buffers.emplace(ev->cpid().get_proc_seq(),
                      std::unordered_map<SymAddrSize, std::list<ZEvent *>>());
    if (isWriteB(ev)) {
      if (!buffers[ev->cpid().get_proc_seq()].count(ev->ml))
        buffers[ev->cpid().get_proc_seq()].emplace(ev->ml, std::list<ZEvent *>());
      buffers[ev->cpid().get_proc_seq()][ev->ml].push_back(ev);
    }
    // Set up proper WriteB <-> WriteM pointers
    if (isWriteM(ev)) {
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
    }
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
  mutated_committed.emplace(readlock->id());
  assert((mutation == ZEventID(true) || graph.hasEvent(mutation)) &&
         graph.hasEvent(readlock));
  const ZEvent * mut_ev = mutation == ZEventID(true) ?
    nullptr : graph.getEvent(mutation);
  assert(!mut_ev || isWriteB(mut_ev) || isUnlock(mut_ev));
  auto remaining_proc = mutated_graph.all_proc_seq();
  for (const auto& tauidx_sch : schedules) {
    const ZEvent * causal = ann_trace.tau[tauidx_sch.first].get();
    assert(isRead(causal) || isLock(causal));
    assert(graph.hasEvent(causal));
    if (!remaining_proc.count(causal->cpid().get_proc_seq()))
      continue;
    if (*causal == *readlock ||
        (!graph.getPo().hasEdge(causal, readlock) &&
         (!mut_ev || !graph.getPo().hasEdge(causal, mut_ev)))) {
      assert(remaining_proc.count(causal->cpid().get_proc_seq()));
      remaining_proc.erase(causal->cpid().get_proc_seq());
      if (remaining_proc.empty())
        break;
    }
    if (!mutated_committed.count(causal->id()))
      mutated_committed.emplace(causal->id());
  }
  time_copy += (double)(clock() - init)/CLOCKS_PER_SEC;
  // Add successful schedule
  assert(schedules.count(readlock->trace_id()) &&
         failed_schedules.count(readlock->trace_id()));
  assert(!schedules.at(readlock->trace_id()).count(mutated_key) &&
         !failed_schedules.at(readlock->trace_id()).count(mutated_key));
  schedules[readlock->trace_id()].emplace(mutated_key,
    ZTrace(&ann_trace, std::move(linear), std::move(mutated_tau),
           std::move(mutated_annotation), std::move(mutated_committed)));
  end_err("mutate-done-added");
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
    err_msg(std::string("recur-on-") + std::to_string(i) +
            "-(idx-is-" + std::to_string(idx) + ")...");
    assert(schedules.count(idx));
    // Move the schedules out of the map
    std::map<ZAnnotation, ZTrace> sch = std::move(schedules.at(idx));
    assert(schedules.at(idx).empty());
    assert(failed_schedules.count(idx));
    failed_schedules[idx].clear();
    // Recur on each schedule
    for (auto& key_trace : sch) {
      if (get_extension(key_trace.second))
        return true;
    }
    // We erase the key idx from (failed)_schedules only *after* the
    // recursive calls, so that we keep the idx as a pointer to the
    // specific read/lock for when iterating over all reads/locks
    assert(schedules.at(idx).empty());
    schedules.erase(idx);
    assert(!schedules.count(idx));
    assert(failed_schedules.at(idx).empty());
    failed_schedules.erase(idx);
    assert(!failed_schedules.count(idx));
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
  if (INFO) ext.dump();
  return extend_and_explore(ann_trace, std::move(ext));
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
  std::unordered_map<SymAddrSize, int> lastUnlock;
  std::unordered_set<SymAddrSize> locked;
  // ThreadID -> ML -> Writes in store queue of thr for ml
  std::unordered_map
    <int, std::unordered_map
     <SymAddrSize, std::list<int>>> storeQueue;
  // Real observation in trace
  std::unordered_map<int, int> realObs;
  // BufferWrite -> pos
  std::unordered_map<ZEventID, unsigned> bw_pos;
  // BufferWrite -> MemoryWrite
  std::unordered_map<int, int> buf_mem;

  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = &(trace.at(i));
    ev->_thread_id = parentTrace.graph.getThreadIDnoAdd(ev);
    assert(ev->thread_id() < 1337);

    if (isWriteB(ev)) {
      bw_pos.emplace(ev->id(), i);
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
      // Check whether the lock is consistent with the annotation
      assert(annotation.lock_defines(ev));
      const auto& an = annotation.lock_obs(ev);
      if (!lastUnlock.count(ev->ml))
        assert(an == ZEventID(true)); // initial
      else
        assert(an == trace.at(lastUnlock.at(ev->ml)).id());
    }
    if (isUnlock(ev)) {
      assert(locked.count(ev->ml));
      locked.erase(ev->ml);
      lastUnlock[ev->ml] = i;
    }
  }

  for (unsigned i=0; i<trace.size(); ++i) {
    const ZEvent *ev = &(trace.at(i));
    if (isRead(ev)) {
      assert(annotation.defines(ev->id()));
      const ZEventID& obs = annotation.obs(ev->id());
      const ZEvent *obsB = parentTrace.graph.initial();
      if (obs.event_id() >= 0) {
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

  // Unset thread-id, let Graph take care of it
  for (const ZEvent& evref : trace)
    evref._thread_id = 1337;

  return true;
}

*/
