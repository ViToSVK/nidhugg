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

#include "ZLinNoclosure.h"

static const bool DEBUG = false;
#include "ZDebug.h"


#define KEY_TSO KeyTSO
#define KEY_PSO MainReqsKeyPSO

#ifndef KEY_TSO
  #define KEY_TSO DummyKey
#endif

#ifndef KEY_PSO
  #define KEY_PSO DummyKey
#endif


/* *************************** */
/* COMMON COMPARISON OPERATORS */
/* *************************** */


template<class T>
bool operator< (const std::vector<T>& left, const std::vector<T>& right) {
  unsigned n = left.size();
  if (n != right.size()) {
    return n < right.size();
  }
  for (unsigned i = 0; i < n; i++) {
    if (left.at(i) != right.at(i)) {
      return left.at(i) < right.at(i);
    }
  }
  return false;
}

template<class T>
bool operator> (const std::vector<T>& left, const std::vector<T>& right) {
  return right < left;
}

template<class T>
bool operator!= (const std::vector<T>& left, const std::vector<T>& right) {
  return left < right || left > right;
}

template<class T>
bool operator== (const std::vector<T>& left, const std::vector<T>& right) {
  return !(left != right);
}


template<class T>
bool operator< (const std::set<T>& left, const std::set<T>& right) {
  if (left.size() != right.size()) {
    return left.size() < right.size();
  }
  auto it2 = right.begin();
  for (auto it1 = left.begin(); it1 != left.end(); it1++, it2++) {
    if (*it1 != *it2) {
      return *it1 < *it2;
    }
  }
  return false;
}

template<class T>
bool operator> (const std::set<T>& left, const std::set<T>& right) {
  return right < left;
}

template<class T>
bool operator!= (const std::set<T>& left, const std::set<T>& right) {
  return left < right || left > right;
}

template<class T>
bool operator== (const std::set<T>& left, const std::set<T>& right) {
  return !(left != right);
}


template<class K, class V>
bool operator< (const std::map<K, V>& left, const std::map<K, V>& right) {
  unsigned n = left.size();
  if (n != right.size()) {
    return n < right.size();
  }
  auto it2 = right.begin();
  for (auto it1 = left.begin(); it1 != left.end(); it1++, it2++) {
    if (*it1 != *it2) {
      return *it1 < *it2;
    }
  }
  return false;
}

template<class K, class V>
bool operator> (const std::map<K, V>& left, const std::map<K, V>& right) {
  return right < left;
}

template<class K, class V>
bool operator!= (const std::map<K, V>& left, const std::map<K, V>& right) {
  return left < right || left > right;
}

template<class K, class V>
bool operator== (const std::map<K, V>& left, const std::map<K, V>& right) {
  return !(left != right);
}


/* *************************** */
/* GENERAL (both TSO and PSO)  */
/* *************************** */


void ZLinNoclosure::calculateWrMapping() {
  start_err("Calculating wr_mapping...");
  for (auto it = an.begin(); it != an.end(); it++) {
    assert(it->first.event_id() >= 0); // read cannot be initial event
    unsigned thr = gr.getThreadIDnoAdd(it->first.cpid().get_proc_seq());
    assert(thr < 10000 && gr.thr_to_lin_id.count(thr));
    thr = gr.thr_to_lin_id.at(thr);
    assert(thr < gr.number_of_threads());
    unsigned ev = (unsigned) it->first.event_id();
    ZObs obsR(thr, ev);
    //const ZObs& obsR = it->first;
    thr = INT_MAX;
    ev = INT_MAX;
    if (it->second.event_id() >= 0) {
      // not initial
      thr = gr.getThreadIDnoAdd(it->second.cpid().get_proc_seq());
      assert(thr < 10000 && gr.thr_to_lin_id.count(thr));
      thr = gr.thr_to_lin_id.at(thr);
      assert(thr < gr.number_of_threads());
      ev = (unsigned) it->second.event_id();
    }
    ZObs obsW(thr, ev);
    //const ZObs& obsW = it->second;
    if (!obsW.isInitial()) {
      wr_mapping[obsW].insert(obsR);
    }
    else {
      const ZEvent *evR = gr.getEvent(it->first); // getEvent(obsR)
      wr_initial[evR->ml].insert(obsR);
    }
  }
  end_err();
}


const ZLinNoclosure::WrSet& ZLinNoclosure::initialGetObservers(SymAddrSize ml) const {
  start_err(std::string("initialGetObservers_") + ml.to_string() + "...");
  auto it = wr_initial.find(ml);
  if (it == wr_initial.end()) {
    end_err("0");
    return dummy;
  }
  end_err("1");
  return it->second;
}


const ZLinNoclosure::WrSet& ZLinNoclosure::getObservers(const ZObs& obs) const {
  start_err(std::string("getObservers_") + obs.to_string() + "...");
  auto it = wr_mapping.find(obs);
  if (it == wr_mapping.end()) {
    end_err("0");
    return dummy;
  }
  end_err("1");
  return it->second;
}
const ZLinNoclosure::WrSet& ZLinNoclosure::getObservers(const ZEvent *ev) const {
  return getObservers(ZObs(toWriteB(ev), gr));
}


unsigned ZLinNoclosure::numEventsInThread(unsigned thr, int aux) const {
  start_err(std::string("numEventsInThread_thr") + std::to_string(thr) + "_aux" + std::to_string(aux) + "...");
  assert(thr < gr.number_of_threads() && "Non-existent thread");
  if (!gr.hasThreadAux(gr.lin_to_thr_id.at(thr), aux)) {
    end_err("0");
    return 0;
  }
  assert(gr.auxes(gr.lin_to_thr_id.at(thr)).count(aux) && "Non-existent aux for given thread");
  LineT line = gr(gr.lin_to_thr_id.at(thr), aux);
  unsigned res = line.size();
  end_err("1");
  return res;
}


/* *************************** */
/* ZSTATE                      */
/* *************************** */


const ZEvent * ZLinNoclosure::State::currEvent(unsigned thr, int aux) const {
  start_err(std::string("currEvent_thr") + std::to_string(thr) + "_aux" + std::to_string(aux) + "...");
  assert(thr < par.gr.number_of_threads() && "Non-existent main thread");
  unsigned pos = prefix.at(thr, aux);
  if (pos >= par.numEventsInThread(thr, aux)) {
    end_err("0");
    return nullptr;
  }
  auto res = par.gr.getEvent(gr.lin_to_thr_id.at(thr), aux, pos);
  end_err("1");
  return res;
}


void ZLinNoclosure::State::add_into_queues(const ZEvent *ev) {
  assert(isWriteB(ev));
  assert(!ev->cpid().is_auxiliary());
  // pso
  if (!pso_queue.count(ev->cpid())) {
    pso_queue.emplace(ev->cpid(), std::unordered_map<SymAddrSize, std::list<ZObs>>());
  }
  assert(pso_queue.count(ev->cpid()));
  if (!pso_queue.at(ev->cpid()).count(ev->ml)) {
    pso_queue.at(ev->cpid()).emplace(ev->ml, std::list<ZObs>());
  }
  assert(pso_queue.at(ev->cpid()).count(ev->ml));
  pso_queue.at(ev->cpid()).at(ev->ml).push_back(ZObs(ev, gr));
  // tso not needed to be maintained
}


void ZLinNoclosure::State::remove_from_queues(const ZEvent *ev) {
  assert(isWriteM(ev));
  const ZEvent *buffer = ev->write_other_ptr;
  assert(isWriteB(buffer));
  assert(!buffer->cpid().is_auxiliary());
  // pso
  assert(pso_queue.count(buffer->cpid()));
  assert(pso_queue.at(buffer->cpid()).count(buffer->ml));
  assert(pso_queue.at(buffer->cpid()).at(buffer->ml).front() == ZObs(buffer, gr));
  pso_queue.at(buffer->cpid()).at(buffer->ml).pop_front();
  // tso not needed to be maintained
}


ZObs ZLinNoclosure::State::what_would_read_observe(const ZEvent *ev) const {
  assert(isRead(ev));
  assert(!ev->cpid().is_auxiliary());
  // read from buffer
  if (pso_queue.count(ev->cpid()) &&
      pso_queue.at(ev->cpid()).count(ev->ml) &&
      !pso_queue.at(ev->cpid()).at(ev->ml).empty()) {
    return ZObs(pso_queue.at(ev->cpid()).at(ev->ml).back());
  }
  // read from main memory
  assert(!pso_queue.count(ev->cpid()) ||
         !pso_queue.at(ev->cpid()).count(ev->ml) ||
         pso_queue.at(ev->cpid()).at(ev->ml).empty());
  if (curr_vals.count(ev->ml)) {
    return curr_vals.at(ev->ml);
  }
  // read initial event
  ZObs initial(INT_MAX, INT_MAX);
  assert(initial.isInitial());
  return initial;
}


bool ZLinNoclosure::State::read_would_observe_what_it_should(const ZEvent *ev) const {
  assert(isRead(ev));
  ZObs would_obs = what_would_read_observe(ev);
  assert(par.an.defines(ev));
  const ZEventID& should_obs_id = par.an.obs(ev);
  if (should_obs_id.event_id() < 0) {
    // initial event
    assert(should_obs_id.event_id() == -1);
    return would_obs.isInitial();
  }
  assert(should_obs_id.event_id() >= 0);
  // noninitial event
  assert(par.gr.hasEvent(should_obs_id));
  const ZEvent * should_buffer = par.gr.getEvent(should_obs_id);
  assert(isWriteB(should_buffer));
  assert(sameMl(should_buffer, ev));
  ZObs should_obs(should_buffer, gr);
  return (would_obs == should_obs);
}


bool ZLinNoclosure::State::isClosedVar(SymAddrSize ml) const {
  start_err(std::string("isClosedVar_") + ml.to_string() + "...");
  auto it = curr_vals.find(ml);
  const std::set<ZObs>& observers = (it == curr_vals.end() ? par.initialGetObservers(ml) : par.getObservers(it->second)).toSet();
  for (auto& r : observers) {
    if (r.ev >= prefix.at(r.thr, -1)) {
      end_err("1");
      return true;
    }
  }
  end_err("0");
  return false;
}


bool ZLinNoclosure::State::canAdvanceAux(unsigned thr, int aux) const {
  start_err(std::string("canAdvanceAux_thr") + std::to_string(thr) + "_aux" + std::to_string(aux) + "...");
  assert(aux != -1 && "CanAdvanceAux can only be called on aux threads");
  assert(thr < par.gr.number_of_threads() && "Non-existent thread");
  const ZEvent *ev = currEvent(thr, aux);
  if (!ev) {
    end_err("0a");
    return false;
  }
  assert(isWriteM(ev));
  if (isWriteM(ev)) {
    if (ev->write_other_ptr->event_id() >= prefix.at(thr)) {
      end_err("0b");
      return false;
    }
    auto res = !isClosedVar(ev->ml);
    end_err("1?");
    return res;
  }
  end_err("1");
  return true;
}


void ZLinNoclosure::State::advance(unsigned thr, int aux, std::vector<ZEvent>& res) {
  start_err(std::string("advance_thr") + std::to_string(thr) + "_aux" + std::to_string(aux) + "...");
  assert(aux == -1 || canAdvanceAux(thr, aux) && "Trying to advance non-advancable aux");
  const ZEvent *ev = currEvent(thr, aux);
  assert(ev && "Trying to advance with nullptr");
  if (isWriteM(ev)) {
    // Update memory
    auto it = curr_vals.find(ev->ml);
    if (it != curr_vals.end()) {
      curr_vals.erase(ev->ml);
    }
    curr_vals.emplace(ev->ml, ZObs(ev->write_other_ptr, gr));
    remove_from_queues(ev);
  }
  if (isWriteB(ev)) {
    add_into_queues(ev);
  }
  // Below assertion fails on PSO if we do not have rule1 edges
  // remote_obsw --> read (our PSO algo assumes and uses that they are there)
  assert(!isRead(ev) || read_would_observe_what_it_should(ev));
  res.push_back(ZEvent(*ev, res.size()));
  prefix.at(thr, aux)++;
  // Update tr_pos
  while (tr_pos < par.tr.size()) {
    const ZEvent& evRef = par.tr.at(tr_pos);
    if (gr.proc_seq_to_thread_id.count(evRef.cpid().get_proc_seq())) {
      int thr_id = gr.proc_seq_to_thread_id.at(evRef.cpid().get_proc_seq());
      assert(gr.thr_to_lin_id.count(thr_id));
      if (gr.hasThreadAux(thr_id, evRef.aux_id())) {
        unsigned lin_id = gr.thr_to_lin_id.at(thr_id);
        if (evRef.event_id() >= prefix.at(lin_id, evRef.aux_id())) {
          break;
        }
      }
    }
    tr_pos++;
  }
  end_err();
}


bool ZLinNoclosure::State::isUseless(const ZEvent *ev) const {
  start_err("isUseless...");
  if (!isWriteM(ev)) {
    end_err("1a");
    return true;
  }
  const auto& wr_set = par.getObservers(ev);
  unsigned numT = wr_set.numThreads();
  if (numT > 1) {
    end_err("0a");
    return false;
  }
  if (numT == 0) {
    end_err("1b");
    return true;
  }
  unsigned thr = wr_set.getOnlyThread();
  if (gr.thr_to_lin_id.at(ev->write_other_ptr->thread_id()) != thr) {
    end_err("0b");
    return false;
  }
  unsigned req = wr_set[thr].last;
  unsigned pos = prefix.at(thr);
  auto res = (pos > req);
  end_err("1?");
  return res;
}


bool ZLinNoclosure::State::canPushUp(unsigned thr, int aux) const {
  start_err(std::string("canPushUp_thr") + std::to_string(thr) + "_aux" + std::to_string(aux) + "...");
  // if aux event
  if (aux != -1) {
    auto res = canAdvanceAux(thr, aux) && isUseless(currEvent(thr, aux));
    end_err("1?");
    return res;
  }
  // else main event
  const ZEvent *ev = currEvent(thr, aux);
  if (!ev) {
    end_err("0a");
    return false;
  }
  for (unsigned thr2 = 0; thr2 < par.gr.number_of_threads(); thr2++) {
    for (int aux2 : par.gr.auxes(gr.lin_to_thr_id.at(thr2))) {
      if (thr == thr2 && aux2 == -1) {
        continue;
      }
      int req = par.po.pred(ev, gr.lin_to_thr_id.at(thr2), aux2).second;
      int pos = prefix.at(thr2, aux2);
      if (req >= pos) {
        end_err("0b");
        return false;
      }
    }
  }
  if (isRead(ev) && !read_would_observe_what_it_should(ev)) {
    return false;
  }
  end_err("1");
  return true;
}


bool ZLinNoclosure::State::allPushedUp() const {
  start_err("allPushedUp...");
  for (unsigned thr = 0; thr < par.gr.number_of_threads(); thr++) {
    for (int aux : par.gr.auxes(gr.lin_to_thr_id.at(thr))) {
      if (canPushUp(thr, aux)) {
        end_err("0");
        return false;
      }
    }
  }
  end_err("1");
  return true;
}


void ZLinNoclosure::State::pushUp(std::vector<ZEvent>& res) {
  start_err("pushUp...");
  bool done = false;
  while (!done) {
    done = true;
    for (unsigned thr = 0; thr < par.gr.number_of_threads(); thr++) {
      for (int aux : par.gr.auxes(gr.lin_to_thr_id.at(thr))) {
        while (canPushUp(thr, aux)) {
          advance(thr, aux, res);
          done = false;
        }
      }
    }
  }
  end_err();
}


bool ZLinNoclosure::State::finished() const {
  start_err("finished...");
  for (unsigned thr = 0; thr < par.gr.number_of_threads(); thr++) {
    unsigned pos = prefix.at(thr);
    unsigned tgt = par.numEventsInThread(thr);
    assert(pos <= tgt);
    if (pos != tgt) {
      end_err("0");
      return false;
    }
  }
  end_err("1");
  return true;
}


void ZLinNoclosure::State::finishOff(std::vector<ZEvent>& res) const {
  for (unsigned thr = 0; thr < par.gr.number_of_threads(); thr++) {
    for (int aux : par.gr.auxes(gr.lin_to_thr_id.at(thr))) {
      unsigned tgt = par.numEventsInThread(thr, aux);
      for (unsigned i = prefix.at(thr, aux); i < tgt; i++) {
        const ZEvent *ev = par.gr(gr.lin_to_thr_id.at(thr), aux).at(i);
        assert(ev && "Line contains a nullptr event");
        res.push_back(ZEvent(*ev, res.size()));
      }
    }
  }
}


/* *************************** */
/* LINEARIZE TSO               */
/* *************************** */


bool ZLinNoclosure::KeyTSO::operator< (const KeyTSO& other) const {
  assert(size() == other.size() && "Can compare only two TSOKeys with same size");
  return vals < other.vals;
}


unsigned ZLinNoclosure::trHintTSO(const State& state) const {
  start_err("trHintTSO...");
  if (tr.empty()) {
    end_err("0");
    return 0;
  }
  if (state.tr_pos == tr.size()) {
    end_err("0a");
    return 0;
  }
  const ZEvent& evRef = tr.at(state.tr_pos);
  if (!isWriteM(&evRef)) {
    end_err("0b");
    return 0;
  }
  assert(gr.proc_seq_to_thread_id.count(evRef.cpid().get_proc_seq()));
  int thr_id = gr.proc_seq_to_thread_id.at(evRef.cpid().get_proc_seq());
  auto res = gr.thr_to_lin_id.at(thr_id);
  end_err("1");
  return res;
}


template<class T>
bool ZLinNoclosure::linearizeTSO(State& curr, std::set<T>& marked, std::vector<ZEvent>& res) const {
  start_err("linearizeTSO/3...");
  if (((double)(std::clock() - start_time)/CLOCKS_PER_SEC)
      > time_limit) {
    exceeded_limit = true;
    res.clear();
    end_err("0to");
    return false;
  }

  // Push-up as much as possible (the boring stuff), then update marked
  // and check for victory
  unsigned orig_size = res.size();
  curr.pushUp(res);
  unsigned pushed_size = res.size();
  err_msg("prefix: " + curr.prefix.str());
  T key(curr);
  if (marked.count(key)) {
    end_err("0a");
    return false;
  }
  marked.insert(key);
  if (curr.finished()) {
    curr.finishOff(res);
    assert(res.size() == gr.events_size());
    end_err("1a");
    return true;
  }
  num_parents++;

  // Now we have choices to make (advance which aux?); try them out
  unsigned n = gr.number_of_threads();
  unsigned start_thr = trHintTSO(curr);
  for (unsigned d = 0; d < n; d++) {
    unsigned thr = (start_thr + d) % n;
    if (!curr.canAdvanceAux(thr)) {
      continue;
    }
    num_children++;
    State next(curr);
    next.advance(thr, 0, res);
    if (linearizeTSO(next, marked, res)) {
      end_err("1b");
      return true;
    }
    if (exceeded_limit) {
      assert(res.empty());
      end_err("0toa");
      return false;
    }
    while (res.size() > pushed_size) {
      res.pop_back();
    }
  }
  while (res.size() > orig_size) {
    res.pop_back();
  }
  end_err("0b");
  return false;
}


template<class T>
std::vector<ZEvent> ZLinNoclosure::linearizeTSO() const
{
  start_time = std::clock();
  start_err("linearizeTSO/0...");
  // po.dump();
  assert(gr.size() > 0);
  State start(*this, gr.number_of_threads());
  std::set<T> marked;
  std::vector<ZEvent> res;
  linearizeTSO<T>(start, marked, res);
  end_err();
  // dump_trace(res);
  elapsed_time = (double)(std::clock() - start_time)/CLOCKS_PER_SEC;
  assert(!exceeded_limit || res.empty());
  assert(res.empty() || res.size() == gr.events_size());
  return res;
}

std::vector<ZEvent> ZLinNoclosure::linearizeTSO() const {
  return linearizeTSO<KEY_TSO>();
}


/* *************************** */
/* LINEARIZE PSO               */
/* *************************** */


ZLinNoclosure::RdyAuxesKeyPSO::RdyAuxesKeyPSO(const State& state)
  : main_prefix(state.prefix.numThreads()), gr(state.gr)
{
  for (unsigned thr = 0; thr < numThreads(); thr++) {
    main_prefix.at(thr) = state.prefix.at(thr);
    for (int aux : state.par.gr.auxes(gr.lin_to_thr_id.at(thr))) {
      const ZEvent *ev = state.currEvent(thr, aux);
      if (ev && !state.isUseless(ev)) {
        ready_auxes.emplace(thr, aux);
      }
    }
  }
}

bool ZLinNoclosure::RdyAuxesKeyPSO::operator< (const RdyAuxesKeyPSO& other) const {
  assert(numThreads() == other.numThreads() && "Can compare only KeyPSOs with same number of threads");
  if (main_prefix != other.main_prefix) {
    return main_prefix < other.main_prefix;
  }
  return ready_auxes < other.ready_auxes;
}


ZLinNoclosure::MainReqsKeyPSO::MainReqsKeyPSO(const State& state)
  : main_prefix(state.prefix.numThreads()), gr(state.gr)
{
  for (unsigned thr = 0; thr < numThreads(); thr++) {
    main_prefix.at(thr) = state.prefix.at(thr);
  }
  for (unsigned thr = 0; thr < numThreads(); thr++) {
    std::map<unsigned, unsigned> reqs;
    // Find the next fence and bwrite (between curr_pos and the fence)
    const ZEvent* next_fence = nullptr;
    bool has_bwrite = false;
    for (unsigned pos = main_prefix.at(thr); pos < state.par.numEventsInThread(thr); pos++) {
      const ZEvent* ev = state.par.gr.getEvent(gr.lin_to_thr_id.at(thr), -1, pos);
      if (ev->fence) {
        next_fence = ev;
        break;
      }
      if (isWriteB(ev)) {
        has_bwrite = true;
      }
    }
    /* If no next fence, leave empty. If has_bwrite, put \bot. Otherwise
     * go through all future mwrite Wm's ancestors of the fence, get their
     * observers, and observers of whatever write is sitting on Wm's variable.
     * Merge. */
    if (!next_fence) {
      continue;
    }
    if (has_bwrite) {
      reqs.emplace(thr, UINT_MAX);
    }
    else {
      for (int aux : state.par.gr.auxes(gr.lin_to_thr_id.at(thr))) {
        if (aux == -1) {
          continue;
        }
        std::pair<const ZEvent *, int> p = state.par.po.pred(next_fence, gr.lin_to_thr_id.at(thr), aux);
        const ZEvent *wm = p.first;
        int last_pred = p.second;
        int curr_pos = state.prefix.at(thr, aux);
        if (last_pred < curr_pos) {
          continue;
        }
        // Collect all the observers
        std::vector<ZObs> observers;
        { // the blocking read
          auto it = state.curr_vals.find(wm->ml);
          const WrSet& wr_set = (it != state.curr_vals.end() ?
            state.par.getObservers(it->second) : state.par.initialGetObservers(wm->ml)
          );
          for (ZObs obs : wr_set.toSet()) {
            observers.push_back(obs);
          }
        }
        { // the future mwrite ancestors of Wm
          for (int pos = last_pred - 1; pos >= curr_pos; pos--) {
            const ZEvent *ev = state.par.gr.getEvent(gr.lin_to_thr_id.at(thr), aux, pos);
            if (!isWriteM(ev)) {
              continue;
            }
            const WrSet& wr_set = state.par.getObservers(ev);
            for (ZObs obs : wr_set.toSet()) {
              observers.push_back(obs);
            }
          }
          // in other threads
          for (unsigned thr2 = 0; thr2 < numThreads(); thr2++) {
            if (thr2 == thr) {
              continue;
            }
            int aux2 = state.par.gr.auxForMl(wm->ml, gr.lin_to_thr_id.at(thr2));
            int last_pred2 = state.par.po.pred(wm, gr.lin_to_thr_id.at(thr2), aux2).second;
            int curr_pos2 = state.prefix.at(thr2, aux2);
            for (int pos2 = last_pred2; pos2 >= curr_pos2; pos2--) {
              const ZEvent *ev = state.par.gr.getEvent(gr.lin_to_thr_id.at(thr2), aux2, pos2);
              if (!isWriteM(ev)) {
                continue;
              }
              const WrSet& wr_set = state.par.getObservers(ev);
              for (ZObs obs : wr_set.toSet()) {
                observers.push_back(obs);
              }
            }
          }
        }
        // merge the observers into reqs
        for (ZObs obs : observers) {
          if (obs.ev < state.prefix.at(obs.thr)) {
            continue;
          }
          auto it = reqs.find(obs.thr);
          unsigned nval = (it == reqs.end() ? obs.ev : std::max(obs.ev, it->second));
          reqs[obs.thr] = nval;
        }
      }
    }
    main_reqs.emplace(thr, reqs);
  }
}

bool ZLinNoclosure::MainReqsKeyPSO::operator< (const MainReqsKeyPSO& other) const {
  assert(numThreads() == other.numThreads() && "Can compare only KeyPSOs with same number of threads");
  if (main_prefix != other.main_prefix) {
    return main_prefix < other.main_prefix;
  }
  return main_reqs < other.main_reqs;
}


bool ZLinNoclosure::canForce(const State& state, unsigned thr) const {
  start_err(std::string("canForce_thr") + std::to_string(thr) + "...");
  // assert(state.allPushedUp() && "Can force only on an all-pushed-up state");
  const ZEvent *ev = state.currEvent(thr);
  if (!ev) {
    end_err("0: null");
    return false;
  }
  // General case
  std::unordered_set<SymAddrSize> occupied;
  for (unsigned thr2 = 0; thr2 < gr.number_of_threads(); thr2++) {
    for (int aux2 : gr.auxes(gr.lin_to_thr_id.at(thr2))) {
      if (aux2 == -1 && thr == thr2) {
        continue;
      }
      int req = po.pred(ev, gr.lin_to_thr_id.at(thr2), aux2).second;
      // main thread
      if (aux2 == -1) {
        if (req >= (int)state.prefix.at(thr2, aux2)) {
          end_err("0: other (main req >= pos)");
          return false;
        }
        continue;
      }
      // aux thread
      int pos = state.prefix.at(thr2, aux2);
      if (pos > req) {
        continue;
      }
      if (!state.canAdvanceAux(thr2, aux2)) {
        end_err("0: other (!canAdvanceAux)");
        return false;
      }
      for (int i = pos; i <= req; i++) {
        const ZEvent *evAux = gr.getEvent(gr.lin_to_thr_id.at(thr2), aux2, i);
        if (occupied.count(evAux->ml)) {
          end_err("0: other (too many advances)");
          return false;
        }
        occupied.insert(evAux->ml);
      }
    }
  }
  end_err("1");
  return true;
}


void ZLinNoclosure::force(State& state, unsigned thr, std::vector<ZEvent>& res) const {
  start_err(std::string("force_thr") + std::to_string(thr) + "...");
  assert(canForce(state, thr) && "According to .canForce, cannot force");
  const ZEvent *ev = state.currEvent(thr);
  for (unsigned thr2 = 0; thr2 < gr.number_of_threads(); thr2++) {
    for (int aux2 : gr.auxes(gr.lin_to_thr_id.at(thr2))) {
      if (aux2 == -1) {
        continue;
      }
      int req = po.pred(ev, gr.lin_to_thr_id.at(thr2), aux2).second;
      while ((int)state.prefix.at(thr2, aux2) <= req) {
        state.advance(thr2, aux2, res);
      }
    }
  }
  state.advance(thr, -1, res);
  end_err();
}


void ZLinNoclosure::calculateTrNextMain() {
  start_err("calculateTrNextMain...");
  int n = tr.size();
  tr_next_main.clear();
  if (tr.empty()) { return; }
  tr_next_main.resize(n+1, UINT_MAX);
  for (int i = n-1; i >= 0; i--) {
    tr_next_main.at(i) = (
      (tr.at(i).aux_id() == -1 &&
       gr.proc_seq_to_thread_id.count(tr.at(i).cpid().get_proc_seq()))
       ? i : tr_next_main.at(i+1));
  }
  end_err();
}


unsigned ZLinNoclosure::trHintPSO(const State& state) const {
  start_err("trHintPSO...");
  if (tr.empty()) {
    end_err("0");
    return 0;
  }
  assert(!tr.empty());
  assert(state.tr_pos < tr_next_main.size());
  unsigned pos = tr_next_main.at(state.tr_pos);
  if (pos == UINT_MAX) {
    end_err("0");
    return 0;
  }
  assert(pos < tr.size());
  const ZEvent& evRef = tr.at(pos);
  assert(gr.proc_seq_to_thread_id.count(evRef.cpid().get_proc_seq()));
  int thr_id = gr.proc_seq_to_thread_id.at(evRef.cpid().get_proc_seq());
  auto res = gr.thr_to_lin_id.at(thr_id);
  end_err("?");
  return res;
}


template<class T>
bool ZLinNoclosure::linearizePSO(State& curr, std::set<T>& marked, std::vector<ZEvent>& res) const {
  start_err("linearizePSO/3...");
  if (((double)(std::clock() - start_time)/CLOCKS_PER_SEC)
      > time_limit) {
    exceeded_limit = true;
    res.clear();
    end_err("0to");
    return false;
  }

  // Push-up as much as possible (the boring stuff), then update marked
  // and check for victory
  unsigned orig_size = res.size();
  curr.pushUp(res);
  unsigned pushed_size = res.size();
  err_msg("prefix: " + curr.prefix.str());
  T key(curr);
  if (marked.count(key)) {
    end_err("0a");
    return false;
  }
  marked.insert(key);
  if (curr.finished()) {
    curr.finishOff(res);
    assert(res.size() == gr.events_size());
    end_err("1a");
    return true;
  }
  num_parents++;

  // Now we have choices to make (which main?); try them out
  unsigned n = gr.number_of_threads();
  unsigned start_thr = trHintPSO(curr);
  for (unsigned d = 0; d < n; d++) {
    unsigned thr = (start_thr + d) % n;
    if (!canForce(curr, thr)) {
      continue;
    }
    num_children++;
    State next(curr);
    force(next, thr, res);
    if (linearizePSO(next, marked, res)) {
      end_err("1b");
      return true;
    }
    if (exceeded_limit) {
      assert(res.empty());
      end_err("0toa");
      return false;
    }
    while (res.size() > pushed_size) {
      res.pop_back();
    }
  }
  while (res.size() > orig_size) {
    res.pop_back();
  }
  end_err("0b");
  return false;
}


template<class T>
std::vector<ZEvent> ZLinNoclosure::linearizePSO() const
{
  start_time = std::clock();
  start_err("linearizePSO/0...");
  // po.dump();
  assert(gr.size() > 0);
  State start(*this, gr.number_of_threads());
  std::set<T> marked;
  std::vector<ZEvent> res;
  linearizePSO<T>(start, marked, res);
  end_err();
  // dump_trace(res);
  elapsed_time = (double)(std::clock() - start_time)/CLOCKS_PER_SEC;
  assert(!exceeded_limit || res.empty());
  if (!(res.empty() || res.size() == gr.events_size())) {
    dump_trace(res);
  }
  assert(res.empty() || res.size() == gr.events_size());
  return res;
}

std::vector<ZEvent> ZLinNoclosure::linearizePSO() const {
  return linearizePSO<KEY_PSO>();
}
