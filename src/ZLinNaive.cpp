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
/* CONSTRUCTOR                 */
/* *************************** */


ZLinNaive::ZLinNaive
(const ZAnnotation& annotation,
 const ZPartialOrder& partialOrder,
 const std::vector<ZEvent>& trace)
: an(annotation), gr(partialOrder.graph),
  po(partialOrder), tr(trace)
{
  calculateWrMapping();
}


void ZLinNaive::calculateWrMapping()
{
  for (auto it = an.begin(); it != an.end(); it++) {
    assert(it->first.event_id() >= 0); // read cannot be initial event
    const ZEvent * read = gr.getEvent(it->first);
    assert(isRead(read));
    if (it->second.event_id() < 0) {
      // initial-event observation
      assert(it->second.event_id() == -1);
      if (!wr_initial.count(read->ml)) {
        wr_initial.emplace(read->ml, std::map<CPid, const ZEvent *>());
      }
      assert(wr_initial.count(read->ml));
      // erase if there is already a read with same cpid and lower event_id
      if (wr_initial.at(read->ml).count(read->cpid())) {
        assert(isRead(wr_initial.at(read->ml).at(read->cpid())));
        if (wr_initial.at(read->ml).at(read->cpid())->event_id() < read->event_id()) {
          wr_initial.at(read->ml).erase(read->cpid());
        }
      }
      // emplace if empty/erased
      if (!wr_initial.at(read->ml).count(read->cpid())) {
        wr_initial.at(read->ml).emplace(read->cpid(), read);
      }
      assert(wr_initial.at(read->ml).count(read->cpid()));
      continue;
    }
    // write observation
    assert(it->second.event_id() >= 0);
    const ZEvent * obsbuf = gr.getEvent(it->second);
    assert(isWriteB(obsbuf));
    const ZEvent * obsmem = obsbuf->write_other_ptr;
    assert(isWriteM(obsmem) && sameMl(obsmem, obsbuf));
    if (!wr_mapping.count(obsmem)) {
      wr_mapping.emplace(obsmem, std::map<CPid, const ZEvent *>());
    }
    assert(wr_mapping.count(obsmem));
    assert(!wr_mapping.count(obsbuf));
    // erase if there is already a read with same cpid and lower event_id
    if (wr_mapping.at(obsmem).count(read->cpid())) {
      assert(isRead(wr_mapping.at(obsmem).at(read->cpid())));
      if (wr_mapping.at(obsmem).at(read->cpid())->event_id() < read->event_id()) {
        wr_mapping.at(obsmem).erase(read->cpid());
      }
    }
    // emplace if empty/erased
    if (!wr_mapping.at(obsmem).count(read->cpid())) {
      wr_mapping.at(obsmem).emplace(read->cpid(), read);
    }
    assert(wr_mapping.at(obsmem).count(read->cpid()));
  }
}


/* *************************** */
/* ZSTATE                      */
/* *************************** */


ZLinNaive::State::State(const ZLinNaive& par0)
 : par(par0), gr(par.gr)
{
  assert(gr.number_of_threads() == gr.get_threads().size());
  assert(gr.number_of_threads() >= 2);
  for (unsigned thr_id : gr.get_threads()) {
    assert(!positions.count(thr_id));
    positions.emplace(thr_id, std::map<int, int>());
    assert(!gr.auxes(thr_id).empty());
    for (int aux_id : gr.auxes(thr_id)) {
      assert(gr.hasThreadAux(thr_id, aux_id));
      // Initialize positions
      assert(!positions.at(thr_id).count(aux_id));
      positions.at(thr_id).emplace(aux_id, 0);
      // Add to next_events_ready(_unique) / next_events_req
      const ZEvent * next = gr.getEvent(thr_id, aux_id, 0);
      assert(next);
      add_to_next_events(next);
    }
  }
}


void ZLinNaive::State::add_to_next_events(const ZEvent * ev)
{
  assert(ev);
  assert(!next_events_ready.count(ev) && !next_events_req.count(ev) &&
         !next_events_ready_unique.count(ev));
  PositionsT req;
  for (unsigned thr_id : gr.get_threads()) {
    for (int aux_id : gr.auxes(thr_id)) {
      if (ev->thread_id() == thr_id && ev->aux_id() == aux_id) {
        continue;
      }
      int pred = par.po.pred(ev, thr_id, aux_id).second;
      if (pred == -1) {
        // No requirement wrt thr_id / aux_id
        continue;
      }
      assert(pred >= 0);
      if (positions.count(thr_id) && positions.at(thr_id).count(aux_id) &&
          positions.at(thr_id).at(aux_id) > pred) {
        // thr_id / aux_id already progressed past the requirement
        continue;
      }
      // Add the thr_id / aux_id requirement of ev
      if (!req.count(thr_id)) {
        req.emplace(thr_id, std::map<int, int>());
      }
      assert(!req.at(thr_id).count(aux_id));
      req.at(thr_id).emplace(aux_id, pred);
    }
  }
  if (req.empty()) {
    if (isRead(ev) || isWriteB(ev) || isWriteM(ev)) {
      next_events_ready.insert(ev);
    } else {
      next_events_ready_unique.insert(ev);
    }
  } else {
    next_events_req.emplace(ev, req);
  }
}


void ZLinNaive::State::update_req(const ZEvent *ev, const std::vector<ZEvent>& res)
{
  assert(ev);
  assert(res.back() == *ev);
  // Update next_events_req
  unsigned thr = ev->thread_id();
  int aux = ev->aux_id();
  int evid = ev->event_id();
  for (auto it = next_events_req.begin(); it != next_events_req.end(); /**/) {
    PositionsT& req = it->second;
    if (!req.count(thr) || !req.at(thr).count(aux) ||
        req.at(thr).at(aux) > evid) {
      ++it;
      continue;
    }
    // Requirement of thr/aux has been met
    assert(req.at(thr).at(aux) == evid);
    req.at(thr).erase(aux);
    if (!req.at(thr).empty()) {
      ++it;
      continue;
    }
    assert(req.at(thr).empty());
    req.erase(thr);
    if (req.empty()) {
      // All requirements have been met, event is ready
      const ZEvent * now_ready = it->first;
      assert(!next_events_ready.count(now_ready) &&
             !next_events_ready_unique.count(now_ready));
      if (isRead(now_ready) || isWriteB(now_ready) || isWriteM(now_ready)) {
        next_events_ready.insert(now_ready);
      } else {
        next_events_ready_unique.insert(now_ready);
      }
      it = next_events_req.erase(it);
      assert(!next_events_req.count(now_ready) &&
             (next_events_ready.count(now_ready) ||
              next_events_ready_unique.count(now_ready)));
    } else {
      ++it;
    }
  }
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


bool ZLinNaive::State::read_was_already_played(const ZEvent *ev) const
{
  assert(isRead(ev));
  assert(positions.count(ev->thread_id()));
  assert(positions.at(ev->thread_id()).count(ev->aux_id()));
  int pos = positions.at(ev->thread_id()).at(ev->aux_id());
  // pos means that events up until including pos-1 have been added
  return (pos > ev->event_id());
}


bool ZLinNaive::State::write_can_overwrite_variable(const ZEvent *ev) const
{
  assert(isWriteM(ev));
  // Make sure ev does not write onto a held variable
  if (main_memory.count(ev->ml)) {
    // Some write has already happened on this ml
    const ZEvent * curr_buff = main_memory.at(ev->ml);
    assert(isWriteB(curr_buff) && sameMl(curr_buff, ev));
    const ZEvent * curr_mem = curr_buff->write_other_ptr;
    assert(isWriteM(curr_mem) && sameMl(curr_mem, curr_buff));
    assert(curr_mem != ev && *curr_mem != *ev);
    assert(!par.wr_mapping.count(curr_buff));
    if (!par.wr_mapping.count(curr_mem)) {
      return true;
    }
    assert(par.wr_mapping.count(curr_mem));
    assert(!par.wr_mapping.at(curr_mem).empty());
    for (const auto& cpid_last : par.wr_mapping.at(curr_mem)) {
      // last read of cpid observing curr_mem
      const ZEvent * last = cpid_last.second;
      assert(isRead(last) && sameMl(last, curr_mem));
      if (!read_was_already_played(last)) {
        // curr_mem cannot be overwritten now, last needs to see it
        return false;
      }
    }
    return true;
  }
  // There is still the initial event on this ml
  assert(!main_memory.count(ev->ml));
  if (!par.wr_initial.count(ev->ml)) {
    return true;
  }
  assert(par.wr_initial.count(ev->ml));
  assert(!par.wr_initial.at(ev->ml).empty());
  for (const auto& cpid_last : par.wr_initial.at(ev->ml)) {
    // last read of cpid observing the initial event
    const ZEvent * last = cpid_last.second;
    assert(isRead(last) && sameMl(last, ev));
    if (!read_was_already_played(last)) {
      // initial event cannot be overwritten now, last needs to see it
      return false;
    }
  }
  return true;
}


bool ZLinNaive::State::can_advance(const ZEvent *ev) const
{
  start_err("canAdvance...");
  assert(ev);
  if (isRead(ev)) {
    end_err("?r");
    return read_would_observe_what_it_should(ev);
  }
  assert(!isRead(ev));
  if (isWriteM(ev)) {
#ifndef NDEBUG
    const ZEvent * buffer = ev->write_other_ptr;
    assert(isWriteB(buffer) && sameMl(buffer, ev));
    assert(buffer->event_id() < positions.at(buffer->thread_id()).at(buffer->aux_id()));
#endif
    end_err("?w");
    return write_can_overwrite_variable(ev);
  }
  assert(!isRead(ev) && !isWriteM(ev));
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
    // Update main_memory and remove_from_queues
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
    // Update add_into_queues
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
  // Update next_events_ready(_unique) and next_events_req
  if (next_events_ready.count(ev)) {
    next_events_ready.erase(ev);
  } else {
    assert(next_events_ready_unique.count(ev));
    next_events_ready_unique.erase(ev);
  }
  update_req(ev, res);
  int pos = positions.at(ev->thread_id()).at(ev->aux_id());
  assert(pos >= 0 && pos <= gr(ev->thread_id(), ev->aux_id()).size());
  if (pos < gr(ev->thread_id(), ev->aux_id()).size()) {
    const ZEvent * next = gr.getEvent(ev->thread_id(), ev->aux_id(), pos);
    assert(next && next->event_id() == pos &&
           !next_events_ready.count(next) && !next_events_req.count(next) &&
           !next_events_ready_unique.count(next));
    add_to_next_events(next);
  }
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
    end_err("0to");
    return false;
  }
  // Check marked
  T key(curr);
  if (marked.count(key)) {
    end_err("0a");
    return false;
  }
  marked.insert(key);
  // Anything unique to extend with?
  if (!curr.next_events_ready_unique.empty()) {
    // There is a non-read-non-write, proceed uniquely
    const ZEvent * ev = *(curr.next_events_ready_unique.begin());
    assert(ev && !isRead(ev) && !isWriteB(ev) && !isWriteM(ev));
    assert(curr.can_advance(ev));
    // Since we proceed uniquely, no need to copy State
    curr.advance(ev, res);
    if (linearize(curr, marked, res)) {
      end_err("1b");
      return true;
    }
    if (exceeded_limit) {
      assert(res.empty());
      end_err("0toa");
      return false;
    }
    // Backtrack, this path does not lead to a witness
    assert(!res.empty() && res.back() == *ev);
    res.pop_back();
    end_err("0b");
    return false;
  }
  // Anything non-unique to extend with?
  assert(curr.next_events_ready_unique.empty());
  if (curr.next_events_ready.empty()) {
    // Successfully finished
    assert(curr.next_events_req.empty());
    assert(res.size() == gr.events_size());
    end_err("1a");
    return true;
  }
  // Order next_events_ready to extend
  int no_of_next = curr.next_events_ready.size();
  assert(no_of_next > 0);
  std::map<int, const ZEvent *> next_ordered;
  for (const ZEvent * ev : curr.next_events_ready) {
    assert(ev);
    assert(isRead(ev) || isWriteB(ev) || isWriteM(ev));
    if (tr.empty()) {
      // order by cpid-hash
      int hsh = ev->cpid().get_hash();
      assert(!next_ordered.count(hsh));
      while (next_ordered.count(hsh)) { hsh++; }
      next_ordered.emplace(hsh, ev);
    } else {
      // order by auxiliary trace
      int key = ev->exec_trace_id;
      assert(key >= 0);
      assert(!next_ordered.count(key));
      while (next_ordered.count(key)) { key++; }
      next_ordered.emplace(key, ev);
    }
  }
  assert(!next_ordered.empty());
  bool has_child = false;
  // Enumerate choices
  for (const auto& id_ev : next_ordered) {
    assert(no_of_next == curr.next_events_ready.size());
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
      end_err("1c");
      return true;
    }
    if (exceeded_limit) {
      assert(res.empty());
      end_err("0tob");
      return false;
    }
    // Backtrack
    assert(!res.empty() && res.back() == *ev);
    res.pop_back();
  }
  // No choice led to a witness
  end_err("0c");
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
