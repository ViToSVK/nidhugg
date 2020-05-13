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

#ifndef __Z_LINEARIZATION_H__
#define __Z_LINEARIZATION_H__

#include <sstream>

#include "ZGraph.h"
#include "ZAnnotation.h"




class ZLinearization {
 private:
  const ZAnnotation& an;
  const ZGraph& gr;
  const ZPartialOrder& po;
  const std::vector<ZEvent>& tr;  // reference-trace for heuristics



  unsigned numEventsInThread(unsigned thr) const;

/*
  class WrEntry {
   public:
    unsigned first, last;

    WrEntry() : first(UINT_MAX), last(0) {}
  };

  class WrSet {
   private:
    std::map<unsigned, WrEntry> thr_map;
    std::set<ZObs> obs_set;

   public:
    const std::set<ZObs>& toSet() const {
      return obs_set;
    }

    void insert(ZObs obs) {
      WrEntry& entry = thr_map[obs.thr];
      if (obs.ev < entry.first) {
        obs_set.erase(ZObs(obs.thr, entry.first));
        entry.first = obs.ev;
        obs_set.insert(obs);
      }
      if (obs.ev > entry.last) {
        obs_set.erase(ZObs(obs.thr, entry.last));
        entry.last = obs.ev;
        obs_set.insert(obs);
      }
    }

    const WrEntry operator[] (unsigned thr) const {
      if (!thr_map.count(thr)) {
        return WrEntry();
      }
      return thr_map.at(thr);
    }

    unsigned numThreads() const {
      return thr_map.size();
    }

    unsigned getOnlyThread() const {
      assert(thr_map.size() == 1 && "Can getOnlyThread only if there is 1 thread");
      return thr_map.begin()->first;
    }
  };

  const WrSet dummy;   // empty set for useless writes
  std::unordered_map<ZObs, WrSet> wr_mapping;
  std::unordered_map<SymAddrSize, WrSet> wr_initial;

  void calculateWrMapping();



  // Returns the observers of (the initial event for variable ml |
  // the event given by (obs | ev)).

  const WrSet& initialGetObservers(SymAddrSize ml) const;
  const WrSet& getObservers(const ZObs& obs) const;
  const WrSet& getObservers(const ZEvent *ev) const;

  class Prefix {
   private:
    std::vector<std::vector<unsigned>> vals;

   public:
    Prefix(unsigned n) : vals(n) {}

    unsigned numThreads() const {
      return vals.size();
    }
    int numAuxes(unsigned thr) const {
      return (thr < vals.size() ? (int)vals.at(thr).size() - 1 : -1);
    }

    unsigned& at (unsigned thr, int aux = -1) {
      std::vector<unsigned>& thr_vals = vals.at(thr);
      unsigned a = aux + 1;
      unsigned nsize = std::max(a + 1, (unsigned)thr_vals.size());
      thr_vals.resize(nsize, 0);
      return thr_vals.at(a);
    }
    const unsigned at (unsigned thr, int aux = -1) const {
      const std::vector<unsigned>& thr_vals = vals.at(thr);
      unsigned a = aux + 1;
      return (a < thr_vals.size() ? thr_vals.at(a) : 0);
    }

    std::string str() const {
      std::stringstream ss;
      for (unsigned thr = 0; thr < numThreads(); thr++) {
        int numA = numAuxes(thr);
        for (int aux = -1; aux < numA; aux++) {
          ss << at(thr, aux) << (aux != numA - 1 ? " " : "");
        }
        ss << " | ";
      }
      return ss.str();
    }
  };
  */

  class State {
   public:
    const ZLinearization& par;
    std::vector<int> key;  //initial size should be number of threads
    std::unordered_map<SymAddrSize,int> occured;
    //each thread last wrie at an ml corresponds to which id
    // size of last_w should be num_threads (have to set)
    std::vector<std::unordered_map<SymAddrSize,ZEventID> > last_w;
    // std::unordered_map<SymAddrSize, ZObs> curr_vals;    // no mapping if initial
    // unsigned tr_pos;

    State(const ZLinearization& par0);

    // Returns the next event in the given thread, or nullptr if there is none.
    const ZEvent * currEvent(unsigned thr) const;

    bool canForce(unsigned thr) const;
  void force(unsigned thr, std::vector<ZEvent>& res) ;

    //bool isClosedVar(SymAddrSize ml) const;
    //bool canAdvance(unsigned thr) const;
    void advance(unsigned thr, std::vector<ZEvent>& res);

    // Heuristic 1
    // bool isUseless(const ZEvent *ev) const;
    // bool canPushUp(unsigned thr, int aux) const;
    // bool allPushedUp() const;
    void pushUp(std::vector<ZEvent>& res);

    // When all main events are already done, what remains is to play
    // the auxiliary events, in (almost) arbitrary order.
     bool finished() const;
    // void finishOff(std::vector<ZEvent>& res) const;
  };
/*
  class DummyKey {
   public:
    DummyKey(const State& state) {}

    bool operator< (const DummyKey& other) const {
      return true;
    }
  };

  // ************* //
  // TSO only      //
  // ************* //

  class KeyTSO {
   private:
    std::vector<unsigned> vals;

    unsigned size() const {
      return vals.size();
    }

   public:
    KeyTSO(State& state) : vals(state.prefix.numThreads()) {
      for (unsigned thr = 0; thr < size(); thr++) {
        vals.at(thr) = state.prefix.at(thr, 0);
      }
    }

    bool operator< (const KeyTSO& other) const;
    bool operator> (const KeyTSO& other) const {
      return other < (*this);
    }
    bool operator== (const KeyTSO& other) const {
      return !((*this) > other) && !((*this) < other);
    }
    bool operator!= (const KeyTSO& other) const {
      return !((*this) == other);
    }
  };

  // Hints the aux thread we should advance.
  unsigned trHintTSO(const State& state) const;

  template<class T>
  bool linearizeTSO(State& curr, std::set<T>& marked, std::vector<ZEvent>& res) const;

  // ************* //
  // PSO only      //
  // ************* //

  class RdyAuxesKeyPSO {
   private:
    std::vector<unsigned> main_prefix;
    std::set<std::pair<unsigned, int>> ready_auxes;

    unsigned numThreads() const {
      return main_prefix.size();
    }

   public:
    RdyAuxesKeyPSO(const State& state);

    bool operator< (const RdyAuxesKeyPSO& other) const;
    bool operator> (const RdyAuxesKeyPSO& other) const {
      return other < *this;
    }
    bool operator== (const RdyAuxesKeyPSO& other) const {
      return !(*this < other) && !(*this > other);
    }
    bool operator!= (const RdyAuxesKeyPSO& other) const {
      return !(*this == other);
    }
  };

  class MainReqsKeyPSO {
   private:
    std::vector<unsigned> main_prefix;
    std::map<unsigned, std::map<unsigned, unsigned>> main_reqs;

    unsigned numThreads() const {
      return main_prefix.size();
    }

   public:
    MainReqsKeyPSO(const State& state);

    bool operator< (const MainReqsKeyPSO& other) const;
    bool operator> (const MainReqsKeyPSO& other) const {
      return other < *this;
    }
    bool operator== (const MainReqsKeyPSO& other) const {
      return !(*this < other) && !(*this > other);
    }
    bool operator!= (const MainReqsKeyPSO& other) const {
      return !(*this == other);
    }
  };

  bool canForce(const State& state, unsigned thr) const;
  void force(State& state, unsigned thr, std::vector<ZEvent>& res) const;

  std::vector<unsigned> tr_next_main;
  void calculateTrNextMain();

  // Hints the next main thread we should force
  unsigned trHintPSO(const State& state) const;

  template<class T>
  bool linearizePSO(State& curr, std::set<T>& marked, std::vector<ZEvent>& res) const;

*/

  /* ********** */
  /* MAIN STUFF */
  /* ********** */

 public:
  ZLinearization(const ZAnnotation& annotation,
                 ZPartialOrder& partialOrder,
                 const std::vector<ZEvent>& trace)
  : an(annotation), gr(partialOrder.graph),
    po(partialOrder), tr(trace)
  {
    //calculateWrMapping();
    //calculateTrNextMain();
  }

  ZLinearization(const ZLinearization& oth) = delete;
  ZLinearization& operator=(ZLinearization& oth) = delete;
  ZLinearization(ZLinearization&& oth) = delete;
  ZLinearization& operator=(ZLinearization&& oth) = delete;

  bool linearize(State& curr, std::set<std::vector<int> >& marked, std::vector<ZEvent>& res) const ;
  std::vector<ZEvent> linearize() const;
/*
  template<class T>
  std::vector<ZEvent> linearizeTSO() const;
  std::vector<ZEvent> linearizeTSO() const;

  template<class T>
  std::vector<ZEvent> linearizePSO() const;
  std::vector<ZEvent> linearizePSO() const;
*/
  // stats
  mutable unsigned num_parents = 0;
  mutable unsigned num_children = 0;
};

#endif // __Z_LINEARIZATION_H__
