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

  std::set<ZObs> dummy;   // empty set for useless writes
  std::unordered_map<ZObs, std::set<ZObs>> wr_mapping;
  std::unordered_map<SymAddrSize, std::set<ZObs>> wr_initial;

  void calculateWrMapping();

  const std::set<ZObs>& initialGetObservers(SymAddrSize ml) const;
  const std::set<ZObs>& getObservers(const ZObs& obs) const;
  const std::set<ZObs>& getObservers(const ZEvent *ev) const;

  unsigned numEventsInThread(unsigned thr, int aux = -1) const;

  class ZPrefix {
   private:
    std::vector<unsigned> vals;

   public:
    ZPrefix(size_t n) : vals(n, 0) {}
    ZPrefix(std::vector<unsigned>& vals0) : vals(vals0) {}

    size_t size() const {
      return vals.size();
    }

    unsigned& at (size_t i) {
      return vals.at(i);
    }
    const unsigned& at (size_t i) const {
      return vals.at(i);
    }

    bool operator< (const ZPrefix& other) const {
      size_t n1 = size(), n2 = other.size();
      size_t n = std::min(n1, n2);
      for (size_t i = 0; i < n; i++) {
        if (vals[i] != other.vals[i]) {
          return vals[i] < other.vals[i];
        }
      }
      return n1 < n2;
    }
    bool operator> (const ZPrefix& other) const {
      return other < *this;
    }
    bool operator== (const ZPrefix& other) const {
      return !(*this < other) && !(*this > other);
    }

    std::string str() const {
      std::stringstream ss;
      if (!vals.empty()) {
        ss << vals[0];
        for (size_t i = 1; i < vals.size(); i++) {
          ss << " " << vals[i];
        }
      }
      return ss.str();
    }
  };

  class ZState {
   public:
    const ZLinearization& par;
    ZPrefix main, aux;
    std::unordered_map<SymAddrSize, ZObs> curr_vals;

    ZState(const ZLinearization& par0, size_t n) : par(par0), main(n), aux(n) {}

    size_t size() const {
      return main.size();
    }

    const ZEvent * currMainEvent(unsigned thr) const;
    const ZEvent * currAuxEvent(unsigned thr) const;

    bool isClosedVar(SymAddrSize ml) const;
    bool canAdvanceAux(unsigned thr) const;
    void advanceAux(unsigned thr, std::vector<ZEvent>& res);

    bool canPushUp(unsigned thr) const;
    void pushUp(std::vector<ZEvent>& res);

    bool finished() const;
    void finishOff(std::vector<ZEvent>& res) const;
  };

  void checkBasis() const;
  bool linearizeTSO(ZState& curr, std::set<ZPrefix>& marked, std::vector<ZEvent>& res) const;

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

  std::vector<ZEvent> linearizeTSO() const;
  std::vector<ZEvent> linearizePSO() const;
};

#endif // __Z_LINEARIZATION_H__
