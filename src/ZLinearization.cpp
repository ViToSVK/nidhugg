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

#include "ZLinearization.h"
static const bool DEBUG = false;
#include "ZDebug.h"


bool operator<(const SymAddrSize &a,const SymAddrSize &b){ 
  return std::make_pair(a.addr,a.size)<std::make_pair(b.addr,b.size);
}

/*




#define KEY_TSO KeyTSO
#define KEY_PSO MainReqsKeyPSO

#ifndef KEY_TSO
  #define KEY_TSO DummyKey
#endif

#ifndef KEY_PSO
  #define KEY_PSO DummyKey
#endif


// *************************** //
// COMMON COMPARISON OPERATORS //
// *************************** //


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


// *************************** //
// GENERAL (both TSO and PSO)  //
// *************************** //


void ZLinearization::calculateWrMapping() {
  start_err("Calculating wr_mapping...");
  for (auto it = an.begin(); it != an.end(); it++) {
    const ZObs& obsR = it->first;
    const ZObs& obsW = it->second;
    if (!obsW.is_initial()) {
      wr_mapping[obsW].insert(obsR);
    }
    else {
      const ZEvent *evR = gr.event(obsR);
      wr_initial[evR->ml].insert(obsR);
    }
  }
  end_err();
}


const ZLinearization::WrSet& ZLinearization::initialGetObservers(SymAddrSize ml) const {
  start_err("initialGetObservers...");
  auto it = wr_initial.find(ml);
  if (it == wr_initial.end()) {
    end_err("0");
    return dummy;
  }
  end_err("1");
  return it->second;
}


const ZLinearization::WrSet& ZLinearization::getObservers(const ZObs& obs) const {
  start_err("getObservers...");
  auto it = wr_mapping.find(obs);
  if (it == wr_mapping.end()) {
    end_err("0");
    return dummy;
  }
  end_err("1");
  return it->second;
}
const ZLinearization::WrSet& ZLinearization::getObservers(const ZEvent *ev) const {
  return getObservers(ZObs(toWriteB(ev)));
}

*/
unsigned ZLinearization::numEventsInThread(unsigned thr) const {
  start_err("numEventsInThread...");
  assert(thr < gr.size() && "Non-existent thread");
  // if (!gr.hasThreadAux(thr, aux)) {
  //   end_err("0");
  //   return 0;
  // }
  //assert(gr.auxes(thr).count(aux) && "Non-existent aux for given thread");
  LineT line = gr(gr.line_id_to_cpid(thr));
  // const ZEvent *lastEv = line.back();
  // bool ahead = (is_read(lastEv) && !an.defines(lastEv)) ||
  //              (is_lock(lastEv) && !an.is_last_lock(lastEv));
  unsigned res = line.size();
  end_err("1");
  return res;
}
/*

// *************************** //
// ZSTATE                      //
// *************************** //

*/
const ZEvent * ZLinearization::State::currEvent(unsigned thr) const {
  start_err("currEvent...");
  assert(thr < par.gr.size() && "Non-existent main thread");
  unsigned pos = key[thr]+1;
  if (pos >= par.numEventsInThread(thr)) {
    end_err("0");
    return nullptr;
  }
  auto cpid= par.gr.line_id_to_cpid(thr);
  auto res = par.gr.event(cpid, pos);
  end_err("1");
  return res;
}
/*

bool ZLinearization::State::isClosedVar(SymAddrSize ml) const {
  start_err("isClosedVar...");
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


bool ZLinearization::State::canAdvance(unsigned thr) const {
  start_err("canAdvance...");
  //assert(aux != -1 && "CanAdvanceAux can only be called on aux threads");
  assert(thr < par.gr.number_of_threads() && "Non-existent thread");
  const ZEvent *ev = currEvent(thr);
  if (!ev) {
    end_err("0a");
    return false;
  }
  if (ev.kind==Kind::WRITE) {
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

*/
ZLinearization::State::State(const ZLinearization& par0): par(par0){
 // last_w.resize(par.gr.size());
  key.resize(par.gr.size(),-1);
   
}
void ZLinearization::State::advance(unsigned thr,  std::vector<ZEvent>& res) {
  start_err("advanceAux...");
 // assert(aux == -1 || canAdvanceAux(thr, aux) && "Trying to advance non-advancable aux");
  const ZEvent *ev = currEvent(thr);
  if (ev->kind==ZEvent::Kind::WRITE) {
    SymAddrSize ml=ev->_ml;
      if(occured[ml]==0){
        occured[ml]=key.size();
        key.push_back(thr);
        // last_w[thr][ml]=ev->_id;
      }
      else{
        key[occured[ml]]=thr;
        // last_w[thr][ml]=ev->_id;
      }
    // Update memory
    // auto it = curr_vals.find(ev->ml);
    // if (it != curr_vals.end()) {
    //   curr_vals.erase(ev->ml);
    // }
    // curr_vals.emplace(ev->ml, ev->write_other_ptr);
    // (key.variables)[ev.ml()] = thr;
  }
  res.push_back(ZEvent(ev));
  key[thr]++;
  // key.positions[thr]++;
  // Update tr_pos
  // while (tr_pos < par.tr.size()) {
  //   const ZEvent& evRef = par.tr.at(tr_pos);
  //   if (evRef.event_id() >= prefix.at(evRef.thread_id(), evRef.aux_id())) {
  //     break;
  //   }
  //   tr_pos++;
  // }
  end_err();
}

/*
bool ZLinearization::State::isUseless(const ZEvent *ev) const {
  start_err("isUseless...");
  if (!is_writeM(ev)) {
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
  if (ev->write_other_ptr->thread_id() != thr) {
    end_err("0b");
    return false;
  }
  unsigned req = wr_set[thr].last;
  unsigned pos = prefix.at(thr);
  auto res = (pos > req);
  end_err("1?");
  return res;
}


bool ZLinearization::State::canPushUp(unsigned thr, int aux) const {
  start_err("canPushUp...");
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
    for (int aux2 : par.gr.auxes(thr2)) {
      if (thr == thr2 && aux2 == -1) {
        continue;
      }
      int req = par.po.pred(ev, thr2, aux2).second;
      int pos = prefix.at(thr2, aux2);
      if (req >= pos) {
        end_err("0b");
        return false;
      }
    }
  }
  end_err("1");
  return true;
}


bool ZLinearization::State::allPushedUp() const {
  start_err("allPushedUp...");
  for (unsigned thr = 0; thr < par.gr.number_of_threads(); thr++) {
    for (int aux : par.gr.auxes(thr)) {
      if (canPushUp(thr, aux)) {
        end_err("0");
        return false;
      }
    }
  }
  end_err("1");
  return true;
}

*/


bool ZLinearization::State::finished() const {
  start_err("finished...");
  for (unsigned thr = 0; thr < par.gr.size(); thr++) {
    unsigned pos = key[thr];
    unsigned tgt = par.numEventsInThread(thr)-1;
    assert(pos <= tgt);
    if (pos != tgt) {
      end_err("0");
      return false;
    }
  }
  end_err("1");
  return true;
}
/*

void ZLinearization::State::finishOff(std::vector<ZEvent>& res) const {
  for (unsigned thr = 0; thr < par.gr.number_of_threads(); thr++) {
    for (int aux : par.gr.auxes(thr)) {
      unsigned tgt = par.numEventsInThread(thr, aux);
      for (unsigned i = prefix.at(thr, aux); i < tgt; i++) {
        const ZEvent *ev = par.gr(thr, aux).at(i);
        res.push_back(ev->copy(res.size(), true));
      }
    }
  }
}


// *************************** //
// LINEARIZE TSO               //
// *************************** //


bool ZLinearization::KeyTSO::operator< (const KeyTSO& other) const {
  assert(size() == other.size() && "Can compare only two TSOKeys with same size");
  return vals < other.vals;
}


unsigned ZLinearization::trHintTSO(const State& state) const {
  start_err("trHintTSO...");
  if (state.tr_pos == tr.size()) {
    end_err("0a");
    return 0;
  }
  const ZEvent& evRef = tr.at(state.tr_pos);
  if (!is_writeM(&evRef)) {
    end_err("0b");
    return 0;
  }
  auto res = evRef.thread_id();
  end_err("1");
  return res;
}


template<class T>
bool ZLinearization::linearizeTSO(State& curr, std::set<T>& marked, std::vector<ZEvent>& res) const {
  start_err("linearizeTSO/3...");

  // Push-up as much as possible (the boring stuff), then update marked
  // and check for victory
  curr.pushUp(res);
  err_msg("prefix: " + curr.prefix.str());
  T key(curr);
  if (marked.count(key)) {
    end_err("0a");
    return false;
  }
  marked.insert(key);
  if (curr.finished()) {
    curr.finishOff(res);
    end_err("1a");
    return true;
  }
  num_parents++;

  // Now we have choices to make (advance which aux?); try them out
  unsigned n = gr.number_of_threads();
  unsigned orig_size = res.size();
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
    while (res.size() > orig_size) {
      res.pop_back();
    }
  }
  end_err("0b");
  return false;
}


template<class T>
std::vector<ZEvent> ZLinearization::linearizeTSO() const
{
  start_err("linearizeTSO/0...");
  // po.dump();
  assert(gr.size() > 0);
  State start(*this, gr.number_of_threads());
  std::set<T> marked;
  std::vector<ZEvent> res;
  linearizeTSO<T>(start, marked, res);
  end_err();
  // dumpTrace(res);
  return res;
}

std::vector<ZEvent> ZLinearization::linearizeTSO() const {
  return linearizeTSO<KEY_TSO>();
}


// *************************** //
// LINEARIZE PSO               //
// *************************** //


ZLinearization::RdyAuxesKeyPSO::RdyAuxesKeyPSO(const State& state)
  : main_prefix(state.prefix.numThreads())
{
  for (unsigned thr = 0; thr < numThreads(); thr++) {
    main_prefix.at(thr) = state.prefix.at(thr);
    for (int aux : state.par.gr.auxes(thr)) {
      const ZEvent *ev = state.currEvent(thr, aux);
      if (ev && !state.isUseless(ev)) {
        ready_auxes.emplace(thr, aux);
      }
    }
  }
}

bool ZLinearization::RdyAuxesKeyPSO::operator< (const RdyAuxesKeyPSO& other) const {
  assert(numThreads() == other.numThreads() && "Can compare only KeyPSOs with same number of threads");
  if (main_prefix != other.main_prefix) {
    return main_prefix < other.main_prefix;
  }
  return ready_auxes < other.ready_auxes;
}


ZLinearization::MainReqsKeyPSO::MainReqsKeyPSO(const State& state)
  : main_prefix(state.prefix.numThreads())
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
      const ZEvent* ev = state.par.gr.event(thr, -1, pos);
      if (ev->fence) {
        next_fence = ev;
        break;
      }
      if (is_writeB(ev)) {
        has_bwrite = true;
      }
    }
    // If no next fence, leave empty. If has_bwrite, put \bot. Otherwise
     * go through all future mwrite Wm's ancestors of the fence, get their
     * observers, and observers of whatever write is sitting on Wm's variable.
     * Merge. //
    if (!next_fence) {
      continue;
    }
    if (has_bwrite) {
      reqs.emplace(thr, UINT_MAX);
    }
    else {
      for (int aux : state.par.gr.auxes(thr)) {
        if (aux == -1) {
          continue;
        }
        std::pair<const ZEvent *, int> p = state.par.po.pred(next_fence, thr, aux);
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
            const ZEvent *ev = state.par.gr.event(thr, aux, pos);
            if (!is_writeM(ev)) {
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
            int aux2 = state.par.gr.auxForMl(wm->ml, thr2);
            int last_pred2 = state.par.po.pred(wm, thr2, aux2).second;
            int curr_pos2 = state.prefix.at(thr2, aux2);
            for (int pos2 = last_pred2; pos2 >= curr_pos2; pos2--) {
              const ZEvent *ev = state.par.gr.event(thr2, aux2, pos2);
              if (!is_writeM(ev)) {
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

bool ZLinearization::MainReqsKeyPSO::operator< (const MainReqsKeyPSO& other) const {
  assert(numThreads() == other.numThreads() && "Can compare only KeyPSOs with same number of threads");
  if (main_prefix != other.main_prefix) {
    return main_prefix < other.main_prefix;
  }
  return main_reqs < other.main_reqs;
}
*/

bool ZLinearization::State::canForce(unsigned thr) const {
  start_err("canForce...");
  // assert(state.allPushedUp() && "Can force only on an all-pushed-up state");
  const ZEvent *ev = currEvent(thr);
  if (!ev) {
    end_err("0: null");
    return false;
  }
  // General case
  // std::unordered_set<SymAddrSize> occupied;
  //check for partial order satisfiability
  for (unsigned thr2 = 0; thr2 < par.gr.size(); thr2++) {
    // for (int aux2 : gr.auxes(thr2)) {
      // if (aux2 == -1 && thr == thr2) {
      //   continue;
      // }
      int req = par.po.pred(ev, par.gr.line_id_to_cpid(thr2)).second;
      // main thread
      // if (aux2 == -1) {
      if (req > (int)key[thr2]) {
        end_err("0: other (main req >= pos)");
        return false;
      }
    }
    //  check for good write satisfiability
    if(ev->kind==ZEvent::Kind::READ){
      SymAddrSize ml=ev->_ml;
      if(occured.at(ml)==0)
        return false;
      unsigned thr_no=key[occured.at(ml)];
      CPid ii=par.gr.line_id_to_cpid(thr_no);
      int evid=par.gr.get_tailw_index( ml, ii, key[thr_no]);
      ZEventID idd=par.gr.event(ii,evid)->_id;
      // ZEventID idd=last_w[thr_no].at(ml);
      if(par.an.mapping.at(ev->_id).goodwrites.find(idd)==par.an.mapping.at(ev->_id).goodwrites.end())
        return false;
    }
    // ZAnn check_set=an.mapping[ev->_id].goodwrites;
    // for(auto eventid : check){
    //   int pos=eventid.event_id();
    //   CPid cpd=pos.cpid();
    //   unsigned thr_pos=gr._cpid_to_line[cpd];

    // }

        // continue;
      // }
      // aux thread
      // int pos = state.prefix.at(thr2, aux2);
      // if (pos > req) {
      //   continue;
      // }
      // if (!state.canAdvanceAux(thr2, aux2)) {
      //   end_err("0: other (!canAdvanceAux)");
      //   return false;
      // }
      // for (int i = pos; i <= req; i++) {
      //   const ZEvent *evAux = gr.event(thr2, aux2, i);
      //   if (occupied.count(evAux->ml)) {
      //     end_err("0: other (too many advances)");
      //     return false;
      //   }
      //   occupied.insert(evAux->ml);
      // }
    // }
  
  end_err("1");
  return true;
}


void ZLinearization::State::force(unsigned thr, std::vector<ZEvent>& res){
  start_err("force...");
  assert(canForce(thr) && "According to .canForce, cannot force");
  const ZEvent *ev = currEvent(thr);
  // for (unsigned thr2 = 0; thr2 < gr.number_of_threads(); thr2++) {
  //   for (int aux2 : gr.auxes(thr2)) {
  //     if (aux2 == -1) {
  //       continue;
  //     }
  //     int req = po.pred(ev, thr2, aux2).second;
  //     while ((int)state.prefix.at(thr2, aux2) <= req) {
  //       state.advance(thr2,res);
  //     }
  //   }
  // }
  advance(thr,res);
  end_err();
}
/*

void ZLinearization::calculateTrNextMain() {
  start_err("calculateTrNextMain...");
  int n = tr.size();
  tr_next_main.clear();
  tr_next_main.resize(n+1, UINT_MAX);
  for (int i = n-1; i >= 0; i--) {
    tr_next_main.at(i) = (tr.at(i).aux_id() == -1 ? i : tr_next_main.at(i+1));
  }
  end_err();
}


unsigned ZLinearization::trHintPSO(const State& state) const {
  start_err("trHintPSO...");
  unsigned pos = tr_next_main.at(state.tr_pos);
  if (pos == UINT_MAX) {
    end_err("0");
    return 0;
  }
  end_err("?");
  return tr.at(pos).thread_id();
}
*/

void ZLinearization::State::pushUp(std::vector<ZEvent>& res) {
  start_err("pushUp...");
  bool done = false;
  while (!done) {
    done = true;
    for (unsigned thr = 0; thr < par.gr.size(); thr++) {
      // if() {
        while (currEvent(thr)->kind==ZEvent::Kind::READ&& canForce(thr)) {
          advance(thr, res);
          done = false;
        }
      
    }
  }
  end_err();
}

//template<class T>
bool ZLinearization::linearize(State& curr, std::set<std::vector<int> >& marked, std::vector<ZEvent>& res) const {
  start_err("linearize/3...");

  // Push-up as much as possible (the boring stuff), then update marked
  // and check for victory
  //Heuristic 1
  curr.pushUp(res);
  err_msg("key: " + curr.key.str());
  //Key key(curr);
  if (marked.count(curr.key)) {
    end_err("0a");
    return false;
  }
  marked.insert(curr.key);
  if (curr.finished()) {
    //curr.finishOff(res);
    end_err("1a");
    return true;
  }
  num_parents++;

  // Now we have choices to make (which main?); try them out
  unsigned n = gr.size();
  unsigned orig_size = res.size();
  //unsigned start_thr = trHintPSO(curr);
  for (unsigned d = 0; d < n; d++) {
    unsigned thr = d;
    if (!curr.canForce(thr)) {
      continue;
    }
    num_children++;
    State next=curr;
    next.force(thr, res);
    if (linearize(next, marked, res)) {
      end_err("1b");
      return true;
    }
    while (res.size() > orig_size) {
      res.pop_back();
    }
  }
  end_err("0b");
  return false;
}


//template<class T>
std::vector<ZEvent> ZLinearization::linearize() const
{
  start_err("linearizePSO/0...");
  // po.dump();
  assert(gr.size() > 0);
  State start(*this);
  std::set<std::vector<int> > marked;
  std::vector<ZEvent> res;
  linearize(start, marked, res);
  end_err();
  // dumpTrace(res);
  return res;
}
/*
std::vector<ZEvent> ZLinearization::linearizePSO() const {
  return linearizePSO<KEY_PSO>();
}


*/
