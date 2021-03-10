/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2020 Viktor Toman
 * Copyright (C) 2020 Truc Lam Bui
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

#include "ZLinNaive.h"

static const bool DEBUG = false;
#include "ZDebug.h"


/* *************************** */
/* KEYNAIVE                    */
/* *************************** */


bool ZLinNaive::KeyNaive::operator< (const KeyNaive& other) const {
  unsigned n = positions.size();
  if (n != other.positions.size()) {
    return n < other.positions.size();
  }
  auto it2 = other.positions.begin();
  for (auto it1 = positions.begin(); it1 != positions.end(); it1++, it2++) {
    if (it1->first != it2->first) {
      return it1->first < it2->first;
    }
    unsigned aux_n = it1->second.size();
    if (aux_n != it2->second.size()) {
      return aux_n < it2->second.size();
    }
    auto aux_it2 = it2->second.begin();
    for (auto aux_it1 = it1->second.begin();
         aux_it1 != it1->second.end(); aux_it1++, aux_it2++) {
      if (*aux_it1 != *aux_it2) {
        return *aux_it1 < *aux_it2;
      }
    }
    assert(aux_it2 == it2->second.end());
  }
  assert(it2 == other.positions.end());
  return false;
}


/* *************************** */
/* ZSTATE                      */
/* *************************** */


ZLinNaive::State::State(const ZLinNaive& par0)
 : par(par0), gr(par.gr)
{
  for (unsigned thr_id : gr.get_threads()) {
    assert(!positions.count(thr_id));
    positions.emplace(thr_id, std::map<int, int>());
    for (int aux_id : gr.auxes(thr_id)) {
      assert(gr.hasThreadAux(thr_id, aux_id));
      assert(!positions.at(thr_id).count(aux_id));
      positions.at(thr_id).emplace(aux_id, 0);
    }
  }
  for (int i = 0; i < par.tr.size(); ++i) {
    int hsh = par.tr.at(i).id().hash();
    assert(!par.event_id_hash_to_trace_id.count(hsh));
    par.event_id_hash_to_trace_id.emplace(hsh, i);
  }
}


// Returns the next events, or empty set if we are finished.
std::unordered_set<const ZEvent *> ZLinNaive::State::next_events() const
{
  start_err("next_events...");
  std::unordered_set<const ZEvent *> res;
  for (unsigned thr_id : gr.get_threads()) {
    assert(positions.count(thr_id));
    for (int aux_id : gr.auxes(thr_id)) {
      assert(gr.hasThreadAux(thr_id, aux_id));
      assert(positions.at(thr_id).count(aux_id));
      int pos = positions.at(thr_id).at(aux_id);
      assert(pos >= 0 && pos <= gr(thr_id, aux_id).size());
      if (pos < gr(thr_id, aux_id).size()) {
        const ZEvent * next = gr.getEvent(thr_id, aux_id, pos);
        assert(next && !res.count(next));
        res.insert(next);
      }
    }
  }
  end_err("1");
  return res;
}


// Returns the PO-minimal events of the input set.
std::unordered_set<const ZEvent *> ZLinNaive::State::po_minimal_events
(const std::unordered_set<const ZEvent *>& input) const
{
  start_err("po_minimal_events...");
  std::unordered_set<const ZEvent *> res(input);
  for (auto it = res.begin(); it != res.end(); ) {
    const ZEvent * ev = *it;
    assert(ev);
    bool minimal = true;
    for (auto it2 = res.begin(); it2 != res.end(); ++it2) {
      if (it == it2) {
        continue;
      }
      const ZEvent * ev2 = *it2;
      assert(ev2);
      assert(ev != ev2 && *ev != *ev2);
      if (par.po.hasEdge(ev2, ev)) {
        minimal = false;
        break;
      }
    }
    if (!minimal) {
      it = res.erase(it);
    } else {
      ++it;
    }
  }
#ifndef NDEBUG
  for (const ZEvent * ev : input) {
    bool minimal = true;
    for (const ZEvent * ev2 : input) {
      if (ev != ev2 && par.po.hasEdge(ev2, ev)) {
        minimal = false;
        break;
      }
    }
    assert(minimal || !res.count(ev));
    assert(!minimal || res.count(ev));
  }
#endif
  end_err("1");
  return res;
}


void ZLinNaive::State::add_into_queues(const ZEvent *ev)
{
  assert(isWriteB(ev));
  assert(!ev->cpid().is_auxiliary());
  // pso
  if (!pso_queue.count(ev->cpid())) {
    pso_queue.emplace(ev->cpid(), std::unordered_map<SymAddrSize, std::list<const ZEvent *>>());
  }
  assert(pso_queue.count(ev->cpid()));
  if (!pso_queue.at(ev->cpid()).count(ev->ml)) {
    pso_queue.at(ev->cpid()).emplace(ev->ml, std::list<const ZEvent *>());
  }
  assert(pso_queue.at(ev->cpid()).count(ev->ml));
  pso_queue.at(ev->cpid()).at(ev->ml).push_back(ev);
  // tso not needed to be maintained
}


void ZLinNaive::State::remove_from_queues(const ZEvent *ev)
{
  assert(isWriteM(ev));
  const ZEvent *buffer = ev->write_other_ptr;
  assert(isWriteB(buffer));
  assert(!buffer->cpid().is_auxiliary());
  // pso
  assert(pso_queue.count(buffer->cpid()));
  assert(pso_queue.at(buffer->cpid()).count(buffer->ml));
  assert(pso_queue.at(buffer->cpid()).at(buffer->ml).front() == buffer);
  pso_queue.at(buffer->cpid()).at(buffer->ml).pop_front();
  // tso not needed to be maintained
}


const ZEvent * ZLinNaive::State::what_would_read_observe(const ZEvent *ev) const
{
  assert(isRead(ev));
  assert(!ev->cpid().is_auxiliary());
  // read from buffer
  if (pso_queue.count(ev->cpid()) &&
      pso_queue.at(ev->cpid()).count(ev->ml) &&
      !pso_queue.at(ev->cpid()).at(ev->ml).empty()) {
    const ZEvent * res = pso_queue.at(ev->cpid()).at(ev->ml).back();
    assert(isWriteB(res));
    return res;
  }
  // read from main memory
  assert(!pso_queue.count(ev->cpid()) ||
         !pso_queue.at(ev->cpid()).count(ev->ml) ||
         pso_queue.at(ev->cpid()).at(ev->ml).empty());
  if (main_memory.count(ev->ml)) {
    const ZEvent * res = main_memory.at(ev->ml);
    assert(isWriteB(res));
    return res;
  }
  // read initial event
  return gr.initial();
}


bool ZLinNaive::State::read_would_observe_what_it_should(const ZEvent *ev) const
{
  assert(isRead(ev));
  const ZEvent * would_obs = what_would_read_observe(ev);
  assert(isInitial(would_obs) || isWriteB(would_obs));
  assert(par.an.defines(ev));
  const ZEventID& should_obs_id = par.an.obs(ev);
  if (should_obs_id.event_id() < 0) {
    // initial event
    assert(should_obs_id.event_id() == -1);
    return isInitial(would_obs);
  }
  assert(should_obs_id.event_id() >= 0);
  // noninitial event
  assert(par.gr.hasEvent(should_obs_id));
  const ZEvent * should_buffer = par.gr.getEvent(should_obs_id);
  assert(isWriteB(should_buffer));
  assert(sameMl(should_buffer, ev));
  return (would_obs == should_buffer);
}


bool ZLinNaive::State::can_advance(const ZEvent *ev) const
{
  start_err("canAdvance...");
  assert(ev);
  if (isRead(ev) && !read_would_observe_what_it_should(ev)) {
    return false;
  }
#ifndef NDEBUG
  if (isWriteM(ev)) {
    const ZEvent * buffer = ev->write_other_ptr;
    assert(isWriteB(buffer) && sameMl(buffer, ev));
    assert(buffer->event_id() < positions.at(buffer->thread_id()).at(buffer->aux_id()));
  }
#endif
  end_err("1");
  return true;
}


void ZLinNaive::State::advance(const ZEvent *ev, std::vector<ZEvent>& res)
{
  assert(ev);
  start_err(std::string("advance_thr") + std::to_string(ev->thread_id())
            + "_aux" + std::to_string(ev->aux_id()) + "...");
  assert(can_advance(ev));
  if (isWriteM(ev)) {
    // Update main memory and queue
    auto it = main_memory.find(ev->ml);
    if (it != main_memory.end()) {
      main_memory.erase(it);
    }
    const ZEvent * buffer = ev->write_other_ptr;
    assert(isWriteB(buffer) && sameMl(buffer, ev));
    main_memory.emplace(ev->ml, buffer);
    remove_from_queues(ev);
  }
  if (isWriteB(ev)) {
    // Update queue
    add_into_queues(ev);
  }
  assert(!isRead(ev) || read_would_observe_what_it_should(ev));
  res.push_back(ZEvent(*ev, res.size()));
  // Update positions
  assert(positions.count(ev->thread_id()));
  assert(positions.at(ev->thread_id()).count(ev->aux_id()));
  assert(positions.at(ev->thread_id()).at(ev->aux_id()) == ev->event_id());
  positions.at(ev->thread_id()).at(ev->aux_id())++;
  assert(positions.at(ev->thread_id()).at(ev->aux_id()) == ev->event_id() + 1);
  end_err();
}


/* *************************** */
/* LINEARIZE                   */
/* *************************** */


template<class T>
bool ZLinNaive::linearize(State& curr, std::set<T>& marked, std::vector<ZEvent>& res) const {
  start_err("linearize/1...");
  if (((double)(std::clock() - start_time)/CLOCKS_PER_SEC)
      > time_limit) {
    exceeded_limit = true;
    res.clear();
    return false;
  }
  // Check marked
  T key(curr);
  if (marked.count(key)) {
    end_err("0a");
    return false;
  }
  marked.insert(key);
  // Anything to extend?
  std::unordered_set<const ZEvent *> next = curr.next_events();
  if (next.empty()) {
    // Successfully finished
    end_err("1a");
    return true;
  }
  // Order partial-order-minimal events to extend
  std::unordered_set<const ZEvent *> next_pomin = curr.po_minimal_events(next);
  assert(!next_pomin.empty());
  std::map<int, const ZEvent *> next_ordered;
  for (const ZEvent * ev : next_pomin) {
    assert(ev);
    if (tr.empty()) {
      assert(!next_ordered.count(ev->cpid().get_hash()));
      next_ordered.emplace(ev->cpid().get_hash(), ev);
    } else {
      int hsh = ev->id().hash();
      assert(event_id_hash_to_trace_id.count(hsh));
      int tr_id = event_id_hash_to_trace_id.at(hsh);
      assert(tr.at(tr_id) == *ev);
      assert(!next_ordered.count(tr_id));
      next_ordered.emplace(tr_id, ev);
    }
  }
  assert(!next_ordered.empty());
  bool has_child = false;
  // Enumerate choices
  for (const auto& id_ev : next_ordered) {
    const ZEvent * ev = id_ev.second;
    if (!curr.can_advance(ev)) {
      continue;
    }
    // Possible choice
    if (!has_child) {
      has_child = true;
      num_parents++;
    }
    num_children++;
    // Fork, advance and recur
    State next(curr);
    next.advance(ev, res);
    if (linearize(next, marked, res)) {
      end_err("1b");
      return true;
    }
    // Backtrack
    assert(res.back() == *ev);
    res.pop_back();
  }
  // No choice led to a witness
  end_err("0b");
  return false;
}


std::vector<ZEvent> ZLinNaive::linearize() const {
  start_time = std::clock();
  start_err("linearize/0...");
  // po.dump();
  assert(!gr.empty() && gr.size() > 0 && gr.number_of_threads() > 1);
  State start(*this);
  std::set<KeyNaive> marked;
  std::vector<ZEvent> res;
  linearize<KeyNaive>(start, marked, res);
  end_err();
  // dump_trace(res);
  elapsed_time = (double)(std::clock() - start_time)/CLOCKS_PER_SEC;
  assert(!exceeded_limit || res.empty());
  return res;
}
