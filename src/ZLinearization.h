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

#ifndef __Z_LINEARIZATION_H__
#define __Z_LINEARIZATION_H__

#include <sstream>

#include "ZGraph.h"


class ZLinearization {
 private:
  const ZAnnotation& an;
  const ZGraph& gr;
  const ZBasis& ba;
  const ZPartialOrder& po;
  const std::vector<ZEvent>& tr;  // reference-trace for heuristics

  unsigned numEventsInThread(unsigned thr, int aux = -1) const;

  std::set<ZObs> dummy;   // empty set for useless writes
  std::unordered_map<ZObs, std::set<ZObs>> wr_mapping;
  std::unordered_map<SymAddrSize, std::set<ZObs>> wr_initial;

  void calculateWrMapping();

  /* Returns the observers of (the initial event for variable ml |
   * the event given by (obs | ev)). */
  const std::set<ZObs>& initialGetObservers(SymAddrSize ml) const;
  const std::set<ZObs>& getObservers(const ZObs& obs) const;
  const std::set<ZObs>& getObservers(const ZEvent *ev) const;

  class ZPrefix {
   private:
    std::vector<std::vector<unsigned>> vals;

   public:
    ZPrefix(unsigned n) : vals(n) {}

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
        for (int aux = -1; aux < numAuxes(thr); aux++) {
          ss << at(thr, aux) << " ";
        }
        ss << "|";
      }
      return ss.str();
    }
  };

  class ZState {
   public:
    const ZLinearization& par;
    ZPrefix prefix;
    std::unordered_map<SymAddrSize, ZObs> curr_vals;
    unsigned tr_pos;

    ZState(const ZLinearization& par0, unsigned n)
      : par(par0), prefix(n), tr_pos(0) {}

    // Returns the next event in the given thread, or nullptr if there is none.
    const ZEvent * currEvent(unsigned thr, int aux = -1) const;

    bool isClosedVar(SymAddrSize ml) const;
    bool canAdvanceAux(unsigned thr, int aux = 0) const;
    void advance(unsigned thr, int aux, std::vector<ZEvent>& res);

    /* What can we play "for free"? */
    bool isUseless(const ZEvent *ev) const;
    bool canPushUp(unsigned thr, int aux) const;
    void pushUp(std::vector<ZEvent>& res);

    /* When all main events are already done, what remains is to play
     * the auxiliary events, in (almost) arbitrary order. */
    bool finished() const;
    void finishOff(std::vector<ZEvent>& res) const;
  };

  /* ************* */
  /* TSO only      */
  /* ************* */

  class ZKeyTSO {
   private:
    std::vector<unsigned> vals;
    
    unsigned size() const {
      return vals.size();
    }
    
   public:
    ZKeyTSO(ZPrefix& prefix) : vals(prefix.numThreads()) {
      for (unsigned thr = 0; thr < size(); thr++) {
        vals.at(thr) = prefix.at(thr, 0);
      }
    }
    
    bool operator< (const ZKeyTSO& other) const {
      assert(size() == other.size() && "Can compare only two TSOKeys with same size");
      for (unsigned thr = 0; thr < size(); thr++) {
        unsigned val1 = vals.at(thr);
        unsigned val2 = other.vals.at(thr);
        if (val1 != val2) {
          return val1 < val2;
        }
      }
      return false;
    }
    bool operator> (const ZKeyTSO& other) const {
      return other < (*this);
    }
    bool operator== (const ZKeyTSO& other) const {
      return !((*this) > other) && !((*this) < other);
    }
    bool operator!= (const ZKeyTSO& other) const {
      return !((*this) == other);
    }
  };
  
  // Hints the main threadID whose aux we should advance.
  unsigned trHintTSO(const ZState& state) const;
  
  bool linearizeTSO(ZState& curr, std::set<ZKeyTSO>& marked, std::vector<ZEvent>& res);

  /* ************* */
  /* PSO only      */
  /* ************* */
  
  // TODO


 public:
  ZLinearization(const ZAnnotation& annotation,
                 ZPartialOrder& partialOrder,
                 const std::vector<ZEvent>& trace)
  : an(annotation), gr(partialOrder.basis.graph),
    ba(partialOrder.basis), po(partialOrder),
    tr(trace)
  {
    calculateWrMapping();
  }

  ZLinearization(const ZLinearization& oth) = delete;
  ZLinearization& operator=(ZLinearization& oth) = delete;
  ZLinearization(ZLinearization&& oth) = delete;
  ZLinearization& operator=(ZLinearization&& oth) = delete;

  std::vector<ZEvent> linearizeTSO();
  std::vector<ZEvent> linearizePSO();
  
  // stats
  unsigned num_parents = 0;
  unsigned num_children = 0;
};

#endif // __Z_LINEARIZATION_H__
