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


std::pair<const ZEvent *, int> ZPartialOrder::succ(const ZEvent *from, unsigned to_thread, int to_aux) const
{
  assert(from && "Null pointer event");
  assert(basis.hasEvent(from));
  assert(from != basis.initial() && "Asking succ for the initial node");
  assert(basis.hasThreadAux(to_thread, to_aux) && "Such thread/aux does not exist");
  assert((from->threadID() != to_thread || from->auxID() != to_aux) && "Same thread/aux");

  unsigned fromLine = basis.lineID(from);
  unsigned toLine = basis.lineID(to_thread, to_aux);
  assert(from->eventID() < _succ.at(fromLine).at(toLine).size());
  int succ_evid = _succ[fromLine][toLine][from->eventID()];
  assert(succ_evid >= 0);

  if ((unsigned) succ_evid >= basis[to_thread, to_aux].size())
    return {nullptr, succ_evid};
  return {basis[to_thread, to_aux][succ_evid], succ_evid};
}


std::pair<const ZEvent *, int> ZPartialOrder::pred(const ZEvent *to, unsigned from_thread, int from_aux) const
{
  assert(to && "Null pointer event");
  assert(basis.hasEvent(to));
  assert(to != basis.initial() && "Asking pred for the initial node");
  assert(basis.hasThreadAux(from_thread, from_aux) && "Such thread/aux does not exist");
  assert((to->threadID() != from_thread || to->auxID() != from_aux) && "Same thread/aux");

  unsigned toLine = basis.lineID(to);
  unsigned fromLine = basis.lineID(from_thread, from_aux);
  assert(to->eventID() < _pred.at(toLine).at(fromLine).size());
  int pred_evid = _pred[toLine][fromLine][to->eventID()];
  assert(pred_evid < (int) basis[from_thread, from_aux].size());

  if (pred_evid < 0)
    return {nullptr, pred_evid};
  return {basis[from_thread, from_aux][pred_evid], pred_evid};
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


void ZPartialOrder::addEdge(const ZEvent *from, const ZEvent *to)
{

}


void ZPartialOrder::addEdgeHelp(unsigned li, unsigned li_evx, unsigned lj, unsigned lj_evx)
{

}


void ZPartialOrder::addLine(const ZEvent * ev)
{

}


void ZPartialOrder::addEvent(const ZEvent * ev)
{

}
