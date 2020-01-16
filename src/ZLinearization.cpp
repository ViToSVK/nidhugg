/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2019 Viktor Toman
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


void ZLinearization::calculateWrMapping() {
  start_err("Calculating wr_mapping...");
  for (auto it = an.begin(); it != an.end(); it++) {
    const ZObs& obsR = it->first;
    const ZObs& obsW = it->second;
    if (!obsW.isInitial()) {
      auto it2 = wr_mapping.emplace(obsW, std::set<ZObs>()).first;
      it2->second.insert(obsR);
    }
    else {
      const ZEvent *evR = ba.getEvent(obsR);
      auto it2 = wr_initial.emplace(evR->ml, std::set<ZObs>()).first;
      it2->second.insert(obsR);
    }
  }
  end_err();
}


const std::set<ZObs>& ZLinearization::initialGetObservers(SymAddrSize ml) const {
  start_err("initialGetObservers...");
  auto it = wr_initial.find(ml);
  if (it == wr_initial.end()) {
    end_err("0");
    return dummy;
  }
  end_err("1");
  return it->second;
}


const std::set<ZObs>& ZLinearization::getObservers(const ZObs& obs) const {
  start_err("getObservers...");
  auto it = wr_mapping.find(obs);
  if (it == wr_mapping.end()) {
    end_err("0");
    return dummy;
  }
  end_err("1");
  return it->second;
}
const std::set<ZObs>& ZLinearization::getObservers(const ZEvent *ev) const {
  return getObservers(ZObs(toWriteB(ev)));
}


unsigned ZLinearization::numEventsInThread(unsigned thr, int aux) const {
  start_err("numEventsInThread...");
  assert(thr < ba.number_of_threads() && "Non-existent thread");
  if (!ba.hasThreadAux(thr, aux)) {
    end_err("0");
    return 0;
  }
  assert(ba.auxes(thr).count(aux) && "Non-existent aux for given thread");
  LineT line = ba(thr, aux);
  const ZEvent *lastEv = line.back();
  bool ahead = (isRead(lastEv) && !an.defines(lastEv)) ||
               (isLock(lastEv) && !an.isLastLock(lastEv));
  unsigned res = line.size() - ahead;
  end_err("1");
  return res;
}


/* *************************** */
/* LINEARIZE TSO               */
/* *************************** */


const ZEvent * ZLinearization::ZState::currMainEvent(unsigned thr) const {
  start_err("currMainEvent...");
  assert(thr < main.size() && "Non-existent thread");
  unsigned pos = main.at(thr);
  unsigned total = par.numEventsInThread(thr);
  assert(pos <= total);
  if (pos == par.numEventsInThread(thr)) {
    end_err("0");
    return nullptr;
  }
  auto res = par.ba.getEvent(thr, -1, pos);
  end_err("1");
  return res;
}


const ZEvent * ZLinearization::ZState::currAuxEvent(unsigned thr) const {
  start_err("currAuxEvent...");
  assert(thr < aux.size() && "Non-existent thread");
  unsigned pos = aux.at(thr);
  if (pos >= par.numEventsInThread(thr, 0)) {
    end_err("0");
    return nullptr;
  }
  auto res = par.ba.getEvent(thr, 0, pos);
  end_err("1");
  return res;
}


bool ZLinearization::ZState::isClosedVar(SymAddrSize ml) const {
  start_err("isClosedVar...");
  auto it = curr_vals.find(ml);
  const std::set<ZObs>& observers = (it == curr_vals.end() ? par.initialGetObservers(ml) : par.getObservers(it->second));
  for (auto& r : observers) {
    if (r.ev >= main.at(r.thr)) {
      end_err("1");
      return true;
    }
  }
  end_err("0");
  return false;
}


bool ZLinearization::ZState::canAdvanceAux(unsigned thr) const {
  start_err("canAdvanceAux...");
  assert(thr < size() && "Non-existent thread");
  const ZEvent *ev = currAuxEvent(thr);
  if (!ev) {
    end_err("0a");
    return false;
  }
  if (ev->write_other_ptr->eventID() >= main.at(thr)) {
    end_err("0b");
    return false;
  }
  auto res = !isClosedVar(ev->ml);
  end_err("1?");
  return res;
}


void ZLinearization::ZState::advanceAux(unsigned thr, std::vector<ZEvent>& res) {
  start_err("advanceAux...");
  assert(canAdvanceAux(thr) && "Trying to advance non-advancable (check first!)");
  const ZEvent *ev = currAuxEvent(thr);
  res.push_back(ev->copy(res.size(), true));
  auto it = curr_vals.find(ev->ml);
  if (it != curr_vals.end()) {
    curr_vals.erase(ev->ml);
  }
  curr_vals.emplace(ev->ml, ev->write_other_ptr);
  aux.at(thr)++;
  end_err();
}


bool ZLinearization::ZState::canPushUp(unsigned thr) const {
  start_err("canPushUp...");
  const ZEvent *ev = currMainEvent(thr);
  if (!ev) {
    end_err("0a");
    return false;
  }
  unsigned n = size();
  for (unsigned other = 0; other < n; other++) {
    if (par.ba.hasThreadAux(other, 0)) {
      int req_aux = par.po.pred(ev, other, 0).second;
      int pos_aux = aux.at(other);
      if (req_aux >= pos_aux) {
        end_err("0b");
        return false;
      }
    }
    if (thr != other) {
      int req_main = par.po.pred(ev, other, -1).second;
      int pos_main = main.at(other);
      if (req_main >= pos_main) {
        end_err("0c");
        return false;
      }
    }
  }
  end_err("1");
  return true;
}


void ZLinearization::ZState::pushUp(std::vector<ZEvent>& res) {
  start_err("pushUp...");
  unsigned n = size();
  bool done = false;
  while (!done) {
    done = true;
    for (unsigned thr = 0; thr < n; thr++) {
      while (canPushUp(thr)) {
        const ZEvent *ev = currMainEvent(thr);
        res.push_back(ev->copy(res.size(), true));
        main.at(thr)++;
        done = false;
      }
    }
  }
  end_err();
}


bool ZLinearization::ZState::finished() const {
  start_err("finished...");
  for (unsigned thr = 0; thr < size(); thr++) {
    unsigned pos = main.at(thr);
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


void ZLinearization::ZState::finishOff(std::vector<ZEvent>& res) const {
  for (unsigned thr = 0; thr < size(); thr++) {
    unsigned tgt = par.numEventsInThread(thr, 0);
    for (unsigned i = aux.at(thr); i < tgt; i++) {
      const ZEvent *ev = par.ba(thr, 0).at(i);
      res.push_back(ev->copy(res.size(), true));
    }
  }
}


bool ZLinearization::linearizeTSO(ZState& curr, std::set<ZPrefix>& marked, std::vector<ZEvent>& res) const {
  start_err("linearizeTSO/3...");
  curr.pushUp(res);
  err_msg("main: " + curr.main.str());
  err_msg("aux:  " + curr.aux.str());
  if (marked.count(curr.aux)) {
    end_err("0a");
    return false;
  }
  marked.insert(curr.aux);
  if (curr.finished()) {
    curr.finishOff(res);
    end_err("1a");
    return true;
  }
  unsigned n = ba.number_of_threads();
  unsigned orig_size = res.size();
  for (unsigned thr = 0; thr < n; thr++) {
    if (!curr.canAdvanceAux(thr)) {
      continue;
    }
    ZState next(curr);
    next.advanceAux(thr, res);
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


std::vector<ZEvent> ZLinearization::linearizeTSO() const
{
  // po.dump();
  assert(ba.size() > 0);
  start_err("linearizeTSO/0...");
  ZState start(*this, ba.number_of_threads());
  std::set<ZPrefix> marked;
  std::vector<ZEvent> res;
  linearizeTSO(start, marked, res);
  end_err();
  // dumpTrace(res);
  return res;
}


/* *************************** */
/* LINEARIZE PSO               */
/* *************************** */

std::vector<ZEvent> ZLinearization::linearizePSO() const
{
  // TODO: implementation for general DC

  return gr.linearizePSO(po, an);
}
