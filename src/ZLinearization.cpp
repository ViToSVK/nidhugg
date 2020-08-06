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
using namespace std;

bool operator<(const SymAddrSize &a,const SymAddrSize &b){
  return std::make_pair(a.addr,a.size)<std::make_pair(b.addr,b.size);
}

std::string pr(std::vector<int> key){
  std::string ans="";
  for(int i = 0; i<key.size(); i++){
    ans+=std::to_string(key[i])+" ";
  }
  return ans;
}

unsigned ZLinearization::numEventsInThread(unsigned thr) const {
   // start_err("numEventsInThread...");
  assert(thr < gr.size() && "Non-existent thread");
  auto cp = gr.line_id_to_cpid(thr);
  LineT line = gr(cp);
  unsigned res = line.size();
  const ZEvent* lastEv = gr.event(cp , res-1);
  bool ahead = (lastEv->kind==ZEvent::Kind::READ && !an.defines(lastEv)) ||
               (lastEv->kind==ZEvent::Kind::M_LOCK && !an.is_last_lock(lastEv));
  res=res-ahead;
  // end_err("1");
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
  int pos = key[thr]+1;
  if (pos >= par.numEventsInThread(thr)) {
    end_err("0");
    return nullptr;
  }
  auto cpid= par.gr.line_id_to_cpid(thr);
  auto res = par.gr.event(cpid, pos);
  end_err("1");
  return res;
}

ZLinearization::State::State(const ZLinearization& par0): par(par0){
  key.resize(par.gr.size(),-1);
  //last_w.resize(par.gr.size());

}
void ZLinearization::State::advance(unsigned thr,  std::vector<ZEvent>& res) {
  // start_err("advanceAux...");
  const ZEvent *ev = currEvent(thr);
  if (ev->kind==ZEvent::Kind::WRITE) {
    // SymAddrSize ml=ev->ml();
      if(occured.find(ev->ml())==occured.end()){
        occured[ev->ml()]=key.size();
        key.push_back(thr);
      }
      else{
        key[occured.at(ev->ml())]=thr;
      }
 //     last_w[thr][ev->ml()]=ev->_id;

  }

  // if((ev->kind==ZEvent::Kind::READ && !par.an.defines(ev)) || ((is_lock(ev) && !par.an.is_last_lock(ev))))
  // {

  // }
  // else{
    res.push_back(*ev);
  // }
  key[thr]++;

  // end_err();
}



bool ZLinearization::State::finished() const {

  start_err("key"+pr(key));
  // start_err("finished...");
  for (unsigned thr = 0; thr < par.gr.size(); thr++) {
    int pos = key[thr];
    int tgt = par.numEventsInThread(thr)-1;
    assert(pos <= tgt);
    if (pos != tgt) {
      end_err("0");
      return false;
    }
  }
  // end_err("1");

  start_err("key"+pr(key));
  return true;
}

bool ZLinearization::State::canForce(unsigned thr) const {
  // start_err("canForce...");
  const ZEvent *ev = currEvent(thr);
  if (!ev) {
    end_err("0: null");
    return false;
  }
  // if(ev->kind==ZEvent::Kind::READ && !par.an.defines(ev))
    // return true;

  //check for partial order satisfiability
  for (unsigned thr2 = 0; thr2 < par.gr.size(); thr2++) {

      if(thr2==thr) continue;
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
      // SymAddrSize ml=ev->_ml;
      if(occured.find(ev->ml())==occured.end()){
        ZEventID idd= par.gr.initial()->_id;
        if(par.an.ann(ev->_id).goodwrites.find(idd)==par.an.ann(ev->_id).goodwrites.end())
        return false;
      else return true;
      }
      unsigned thr_no=key[occured.at(ev->ml())];
      CPid ii=par.gr.line_id_to_cpid(thr_no);
      //int evid=par.gr.get_tailw_index( ev->ml(), ii, key[thr_no]);
      const ZEvent* ev1= par.gr.get_tailw(ev->ml(), ii, key[thr_no]);
      ZEventID idd=ev1->_id;
      if(par.an.ann(ev->_id).goodwrites.find(idd)==par.an.ann(ev->_id).goodwrites.end())
        return false;
      // if(thr==2){
      //   ev->dump();
      //   //par.gr.event(ii,evid)->dump();
      //  // last_w[thr_no].at(ev->ml()).dump();
      //   end_err(std::to_string(key[thr_no]));
      //  end_err(std::to_string(thr_no));

      // }
    }
  // end_err("1");
  return true;
}


void ZLinearization::State::force(unsigned thr, std::vector<ZEvent>& res){
  // start_err("force...");
  assert(canForce(thr) && "According to .canForce, cannot force");
  // const ZEvent *ev = currEvent(thr);

  advance(thr,res);
  end_err();
}


void ZLinearization::State::pushUp(std::vector<ZEvent>& res) {
  start_err("pushUp...");
  bool done = false;
  while (!done) {
    done = true;
    for (unsigned thr = 0; thr < par.gr.size(); thr++) {

        while (currEvent(thr) && currEvent(thr)->kind!=ZEvent::Kind::WRITE && canForce(thr)) {
          advance(thr, res);
          done = false;
        }

    }
  }
  // end_err();
}

bool ZLinearization::linearize(State& curr, std::set<std::vector<int> >& marked, std::vector<ZEvent>& res) const {
  // start_err("linearize/3...");

  //Heuristic 1
  // start_err("key"+pr(curr.key));
  curr.pushUp(res);
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
  bool fl=0;
  vector< pair<int,int> > thr_or;
  // Now we have choices to make (which main?); try them out
  unsigned n = gr.size();
  unsigned orig_size = res.size();
  //unsigned start_thr = trHintPSO(curr);
  for (unsigned d = 0; d < n; d++) {
    unsigned thr = d;
    const ZEvent *ev = curr.currEvent(thr);
    if (!ev) {
      end_err("0: null");
      continue;
    }
    else if(ev->kind!=ZEvent::Kind::WRITE){
      continue;
    }
    if (!curr.canForce(thr)) {
      continue;
    }
    thr_or.push_back(make_pair(ev->trace_id(),thr));
    // num_children++;
    fl=1;

  }
  sort(thr_or.begin(),thr_or.end());
  for(int i=0;i<thr_or.size();i++){
    State next=curr;
    unsigned thr=thr_or[i].second;
    num_children++;
    next.force(thr, res);
    if (linearize(next, marked, res)) {
      // end_err("1b");
      return true;
    }
    while (res.size() > orig_size) {
      res.pop_back();
    }
  }
  if(!fl)
    num_parents--;
  // end_err("0b");
  return false;
}

std::vector<ZEvent> ZLinearization::linearize() const
{
  start_err("linearizePSO/0...");
  // po.dump();
  assert(gr.size() > 0);
  State start(*this);
  std::set<std::vector<int> > marked;
  std::vector<ZEvent> res,res2;
  bool suc=linearize(start, marked, res);
  end_err("finished");

  start_err("key"+pr(start.key));
  // dumpTrace(res);
  if(suc){
    end_err("succeeded");
  return res;
  }
else{
  end_err("linearisation not found");
  return res2;
}
}
