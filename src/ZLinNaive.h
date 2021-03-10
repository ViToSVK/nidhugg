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

#ifndef __Z_LINNAIVE_H__
#define __Z_LINNAIVE_H__

#include <sstream>

#include "ZLinearization.h"


class ZLinNaive {
 private:
  const ZAnnotation& an;
  const ZGraph& gr;
  const ZPartialOrder& po;
  const std::vector<ZEvent>& tr;  // reference-trace for heuristics
  mutable std::map<int, int> event_id_hash_to_trace_id;


  class State {
   public:
    const ZLinNaive& par;
    const ZGraph& gr;
    std::map<unsigned, std::map<int, int>> positions;
    std::map<CPid, std::unordered_map<
      SymAddrSize, std::list<const ZEvent *>>> pso_queue;
    std::unordered_map<SymAddrSize, const ZEvent *> main_memory;

    State(const ZLinNaive& par0);
    State(const State&) = default;
    State(State&&) = default;
    State& operator=(const State&) = delete;
    State& operator=(State&&) = delete;

    // Returns the next events, or empty set if we are finished.
    std::unordered_set<const ZEvent *> next_events() const;
    // Returns the PO-minimal events of the input set.
    std::unordered_set<const ZEvent *> po_minimal_events
    (const std::unordered_set<const ZEvent *>& input) const;

    // Check observation of a read
    void add_into_queues(const ZEvent *ev);
    void remove_from_queues(const ZEvent *ev);
    const ZEvent * what_would_read_observe(const ZEvent *ev) const;
    bool read_would_observe_what_it_should(const ZEvent *ev) const;

    bool can_advance(const ZEvent *ev) const;
    void advance(const ZEvent *ev, std::vector<ZEvent>& res);
  };


  class KeyNaive {
   private:
    std::map<unsigned, std::map<int, int>> positions;

   public:
    KeyNaive(State& state) : positions(state.positions) {}

    bool operator< (const KeyNaive& other) const;
    bool operator> (const KeyNaive& other) const {
      return other < (*this);
    }
    bool operator== (const KeyNaive& other) const {
      return !((*this) > other) && !((*this) < other);
    }
    bool operator!= (const KeyNaive& other) const {
      return !((*this) == other);
    }
  };


  template<class T>
  bool linearize(State& curr, std::set<T>& marked, std::vector<ZEvent>& res) const;

  /* ********** */
  /* MAIN STUFF */
  /* ********** */

  mutable double start_time;
 public:
  ZLinNaive(const ZAnnotation& annotation,
            const ZPartialOrder& partialOrder,
            const std::vector<ZEvent>& trace)
  : an(annotation), gr(partialOrder.graph),
    po(partialOrder), tr(trace) {}

  ZLinNaive(const ZLinNaive& oth) = delete;
  ZLinNaive& operator=(ZLinNaive& oth) = delete;
  ZLinNaive(ZLinNaive&& oth) = delete;
  ZLinNaive& operator=(ZLinNaive&& oth) = delete;

  std::vector<ZEvent> linearize() const;

  // stats
  mutable unsigned num_parents = 0;
  mutable unsigned num_children = 0;
  mutable double elapsed_time;
  double time_limit = 60;
  mutable bool exceeded_limit = false;
};

#endif // __Z_LINNAIVE_H__
