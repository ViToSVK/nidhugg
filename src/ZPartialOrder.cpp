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

#include "ZPartialOrder.h"


std::pair<const ZEvent *, int> ZPartialOrder::succ(const ZEvent *from, unsigned to_line) const
{
  assert(from && "Null pointer event");
  assert(basis.hasEvent(from));
  assert(from != basis.initial() && "Asking succ for the initial node");

  unsigned from_line = basis.lineID(from);
  assert(from_line != to_line);
  assert(from->eventID() < _succ.at(from_line).at(to_line).size());
  int succ_evid = _succ[from_line][to_line][from->eventID()];
  assert(succ_evid >= 0);

  assert(to_line < basis.lines.size());
  if ((unsigned) succ_evid >= basis.lines[to_line].size())
    return {nullptr, succ_evid};
  return {basis.lines[to_line][succ_evid], succ_evid};
}


std::pair<const ZEvent *, int> ZPartialOrder::succ(const ZEvent *from, unsigned to_thread, int to_aux) const
{
  assert(from && "Null pointer event");
  assert(basis.hasEvent(from));
  assert(from != basis.initial() && "Asking succ for the initial node");
  assert(basis.hasThreadAux(to_thread, to_aux) && "Such thread/aux does not exist");
  assert((from->threadID() != to_thread || from->auxID() != to_aux) && "Same thread/aux");

  unsigned toLine = basis.lineID(to_thread, to_aux);
  return succ(from, toLine);
}


std::pair<const ZEvent *, int> ZPartialOrder::pred(const ZEvent *to, unsigned from_line) const
{
  assert(to && "Null pointer event");
  assert(basis.hasEvent(to));
  assert(to != basis.initial() && "Asking pred for the initial node");

  unsigned to_line = basis.lineID(to);
  assert(from_line != to_line);
  assert(to->eventID() < _pred.at(to_line).at(from_line).size());
  int pred_evid = _pred[to_line][from_line][to->eventID()];
  assert(from_line < basis.lines.size());
  assert(pred_evid < (int) basis.lines[from_line].size());

  if (pred_evid < 0)
    return {nullptr, pred_evid};
  return {basis.lines[from_line][pred_evid], pred_evid};
}


std::pair<const ZEvent *, int> ZPartialOrder::pred(const ZEvent *to, unsigned from_thread, int from_aux) const
{
  assert(to && "Null pointer event");
  assert(basis.hasEvent(to));
  assert(to != basis.initial() && "Asking pred for the initial node");
  assert(basis.hasThreadAux(from_thread, from_aux) && "Such thread/aux does not exist");
  assert((to->threadID() != from_thread || to->auxID() != from_aux) && "Same thread/aux");

  unsigned fromLine = basis.lineID(from_thread, from_aux);
  return pred(to, fromLine);
}


bool ZPartialOrder::hasEdge(const ZEvent *from, const ZEvent *to) const
{
  assert(from && to && "Null pointer event");
  assert(from != to);
  assert(basis.hasEvent(from) && basis.hasEvent(to));
  if (from == basis.initial())
    return true; // init HB everything
  if (to == basis.initial())
    return false; // nothing HB init
  if (from->threadID() == to->threadID() &&
      from->auxID() == to->auxID())
    return from->eventID() < to->eventID();

  unsigned fromLine = basis.lineID(from);
  unsigned toLine = basis.lineID(from);
  return _succ[fromLine][toLine]
              [from->eventID()] <= (int) to->eventID();
}


bool ZPartialOrder::areOrdered(const ZEvent *ev1, const ZEvent *ev2) const
{
  return hasEdge(ev1, ev2) || hasEdge(ev2, ev1);
}


// Adds an edge between two unordered nodes
// This method maintains:
// 1) line-pair-wise transitivity
// 2) complete transitivity
void ZPartialOrder::addEdge(const ZEvent *from, const ZEvent *to)
{
  assert(from && to && "Null pointer event");
  assert(!areOrdered(from, to));

  // Maintenance of the complete transitivity
  // Collect nodes from different lines with edges:
  // 1) to   li[li_evx]
  // 2) from lj[lj_evx]

  assert(basis.size() == _succ.size() && basis.size() == _pred.size());
  unsigned li = basis.lineID(from);
  unsigned lj = basis.lineID(to);

  std::set<const ZEvent *> before_from, after_to;
  for (unsigned lk = 0; lk<basis.size(); ++lk) {
    if (lk != li && lk != lj) {
      auto maxbefore = pred(from, lk);
      assert(maxbefore.second < (int) basis.lines[lk].size());
      if (maxbefore.first)
        before_from.emplace(maxbefore.first);
      else assert(maxbefore.second == -1);

      auto minafter = succ(to, lk);
      assert(minafter.second >= 0);
      if (minafter.first)
        after_to.emplace(minafter.first);
      else
        assert(minafter.second >= (int) basis.lines[lk].size());
    }
  }

  addEdgeHelp(from, to);

  // TryAdd edges between each of *before_from* and *to*
  for (const auto& bef : before_from) {
    assert(!hasEdge(to, bef) && "Cycle");
    if (!hasEdge(bef, to))
      addEdgeHelp(bef, to);
  }

  // TryAdd edges between *from* and each of *after_to*
  for (const auto& aft : after_to) {
    assert(!hasEdge(aft, from) && "Cycle");
    if (!hasEdge(from, aft))
      addEdgeHelp(from, aft);
  }

  // TryAdd edges between each of *before_from* and each of *after_to*
  // (if they belong to different threads)
  for (auto& bef : before_from)
    for (auto& aft : after_to)
      if (basis.lineID(bef) != basis.lineID(aft)) {
        assert(!hasEdge(aft, bef) && "Cycle");
        if (!hasEdge(bef, aft))
          addEdgeHelp(bef, aft);
      }
}


// Helper method for addEdge,
// maintains line-pair-wise transitivity
void ZPartialOrder::addEdgeHelp(const ZEvent *from, const ZEvent *to)
{
  unsigned li = basis.lineID(from);
  unsigned lj = basis.lineID(to);
  unsigned li_evx = from->eventID();
  unsigned lj_evx = to->eventID();
  assert( _succ[li][lj][li_evx] > (int) lj_evx && // ! li[li_evx] HB lj[lj_evx]
          _succ[lj][li][lj_evx] > (int) li_evx && // ! lj[lj_evx] HB li[li_evx]
          "Tried to add an edge between ordered nodes");
  assert( _pred[lj][li][lj_evx] < (int) li_evx && // ! li[li_evx] HB lj[lj_evx]
          _pred[li][lj][li_evx] < (int) lj_evx && // ! lj[lj_evx] HB li[li_evx]
          "Inconsistent succ/pred vector clocks");

  _succ[li][lj][li_evx] = (int) lj_evx;
  _pred[lj][li][lj_evx] = (int) li_evx;

  // Maintenance of thread-pair-wise transitivity of _succ[li][lj]
  // Everything in li before li_evx also happens-before lj[lj_evx]

  for (int li_before_evx = li_evx - 1;
       li_before_evx >= 0;
       --li_before_evx) {
    if (_succ[li][lj][li_before_evx] <= (int) lj_evx)
      break; // since for smaller indices also <= lj_evx
    _succ[li][lj][li_before_evx] = (int) lj_evx;
  }

  // Maintenance of thread-pair-wise transitivity of _pred[lj][li]
  // Everything in lj after lj_evx also happens-after li[li_evx]

  for (unsigned lj_after_evx = lj_evx + 1;
       lj_after_evx < basis.lines[lj].size();
       ++lj_after_evx) {
    if (_pred[lj][li][lj_after_evx] >= (int) li_evx)
      break; // since for bigger indices also >= li_evx
    _pred[lj][li][lj_after_evx] = (int) li_evx;
  }
}


void ZPartialOrder::addLine(const ZEvent * ev)
{
  // This should be called right after line is added to basis
  assert(ev && "Null pointer event");
  assert(basis.hasEvent(ev));
  assert(ev->eventID() == 0);
  assert(basis.size() == _succ.size() -1 &&
         basis.size() == _pred.size() -1);
  assert(basis.lineID(ev) == _succ.size());

  // Clocks from original lines to new one
  for (unsigned li=0; li<_succ.size(); li++) {
    assert(_succ.size() == _succ[li].size());
    assert(_pred.size() == _pred[li].size());
    _succ[li].push_back(std::vector<int>(basis.lines[li].size(), INT_MAX));
    _pred[li].push_back(std::vector<int>(basis.lines[li].size(), -1));
    assert(basis.size() == _succ[li].size());
    assert(basis.size() == _pred[li].size());
  }

  // Clocks for new line
  unsigned lnew = _succ.size();
  _succ.push_back(std::vector<std::vector<int>>());
  _succ[lnew].reserve(basis.size());
  _pred.push_back(std::vector<std::vector<int>>());
  _pred[lnew].reserve(basis.size());
  assert(basis.size() == _succ.size() && basis.size() == _pred.size());
  assert(basis.lines.at(lnew).size() == 0);
  for (unsigned lold=0; lold<basis.size(); lold++) {
    _succ[lnew].push_back(std::vector<int>());
    _pred[lnew].push_back(std::vector<int>());
  }
}


void ZPartialOrder::addEvent(const ZEvent * ev)
{
  assert(ev && "Null pointer event");
  assert(basis.hasEvent(ev));
  unsigned lID = basis.lineID(ev);
  unsigned evID = ev->eventID();

  // This should be called right after ev is added to basis
  assert(basis.size() == _succ.size() && basis.size() == _pred.size());
  assert(evID == basis.lines[lID].size() - 1);
  for (unsigned li=0; li<basis.size(); li++) {
    assert(evID == _succ[lID][li].size());
    assert(evID == _pred[lID][li].size());
    // New _succ slot should have INT_MAX,
    _succ[lID][li].push_back(INT_MAX);
    // New _pred slot should have what the last slot says (if there is any; else -1)
    int newpred = (evID > 0) ? _pred[lID][li][evID-1] : -1;
    _pred[lID][li].push_back(newpred);
    assert(basis.lines[lID].size() == _succ[lID][li].size());
    assert(basis.lines[lID].size() == _pred[lID][li].size());
  }
}
