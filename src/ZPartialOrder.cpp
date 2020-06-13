/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2020 Viktor Toman
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
#include "ZGraph.h"


ZGraph ZPartialOrder::graphDummy = ZGraph();


// Initial
ZPartialOrder::ZPartialOrder(const ZGraph& graph)
  : graph(graph),
    _succ(),
    _pred()
{
  assert(empty());
}


// Empty
ZPartialOrder::ZPartialOrder()
  : ZPartialOrder(ZPartialOrder::graphDummy) {}


std::pair<const ZEvent *, int> ZPartialOrder::succ(const ZEvent *from, unsigned to_line) const
{
  assert(from && "Null pointer event");
  assert(graph.hasEvent(from));
  assert(from != graph.initial() && "Asking succ for the initial node");

  unsigned from_line = graph.lineID(from);
  assert(from_line != to_line);
  assert(from->event_id() < _succ.at(from_line).at(to_line).size());
  int succ_evid = _succ[from_line][to_line][from->event_id()];
  assert(succ_evid >= 0);

  assert(to_line < graph.lines.size());
  if ((unsigned) succ_evid >= graph.lines[to_line].size()) {
    assert(succ_evid == INT_MAX);
    return {nullptr, succ_evid};
  }
  return {graph.lines[to_line][succ_evid], succ_evid};
}


std::pair<const ZEvent *, int> ZPartialOrder::succ(const ZEvent *from, unsigned to_thread, int to_aux) const
{
  assert(from && "Null pointer event");
  assert(graph.hasEvent(from));
  assert(from != graph.initial() && "Asking succ for the initial node");
  assert(graph.hasThreadAux(to_thread, to_aux) && "Such thread/aux does not exist");
  assert((from->thread_id() != to_thread || from->aux_id() != to_aux) && "Same thread/aux");

  unsigned toLine = graph.lineID(to_thread, to_aux);
  return succ(from, toLine);
}


std::pair<const ZEvent *, int> ZPartialOrder::pred(const ZEvent *to, unsigned from_line) const
{
  assert(to && "Null pointer event");
  assert(graph.hasEvent(to));
  assert(to != graph.initial() && "Asking pred for the initial node");

  unsigned to_line = graph.lineID(to);
  assert(from_line != to_line);
  assert(to->event_id() < _pred.at(to_line).at(from_line).size());
  int pred_evid = _pred[to_line][from_line][to->event_id()];
  assert(from_line < graph.lines.size());
  assert(pred_evid < (int) graph.lines[from_line].size());

  if (pred_evid < 0) {
    assert(pred_evid == -1);
    return {nullptr, pred_evid};
  }
  return {graph.lines[from_line][pred_evid], pred_evid};
}


std::pair<const ZEvent *, int> ZPartialOrder::pred(const ZEvent *to, unsigned from_thread, int from_aux) const
{
  assert(to && "Null pointer event");
  assert(graph.hasEvent(to));
  assert(to != graph.initial() && "Asking pred for the initial node");
  assert(graph.hasThreadAux(from_thread, from_aux) && "Such thread/aux does not exist");
  assert((to->thread_id() != from_thread || to->aux_id() != from_aux) && "Same thread/aux");

  unsigned fromLine = graph.lineID(from_thread, from_aux);
  return pred(to, fromLine);
}


bool ZPartialOrder::hasEdge(const ZEvent *from, const ZEvent *to) const
{
  assert(from && to && "Null pointer event");
  assert(from != to);
  assert(graph.hasEvent(from) && graph.hasEvent(to));
  if (from == graph.initial())
    return true; // init HB everything
  if (to == graph.initial())
    return false; // nothing HB init
  if (from->thread_id() == to->thread_id() &&
      from->aux_id() == to->aux_id())
    return from->event_id() < to->event_id();

  unsigned fromLine = graph.lineID(from);
  unsigned toLine = graph.lineID(to);
  return _succ[fromLine][toLine]
              [from->event_id()] <= (int) to->event_id();
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

  assert(graph.size() == _succ.size() && graph.size() == _pred.size());
  unsigned li = graph.lineID(from);
  unsigned lj = graph.lineID(to);
  assert(li != lj);

  std::set<const ZEvent *> before_from, after_to;
  for (unsigned lk = 0; lk<graph.size(); ++lk) {
    if (lk != li && lk != lj) {
      auto maxbefore = pred(from, lk);
      assert(maxbefore.second < (int) graph.lines[lk].size());
      if (maxbefore.first)
        before_from.emplace(maxbefore.first);
      else assert(maxbefore.second == -1);

      auto minafter = succ(to, lk);
      assert(minafter.second >= 0);
      if (minafter.first)
        after_to.emplace(minafter.first);
      else
        assert(minafter.second >= (int) graph.lines[lk].size());
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
      if (graph.lineID(bef) != graph.lineID(aft)) {
        assert(!hasEdge(aft, bef) && "Cycle");
        if (!hasEdge(bef, aft))
          addEdgeHelp(bef, aft);
      }
}


// Helper method for addEdge,
// maintains line-pair-wise transitivity
void ZPartialOrder::addEdgeHelp(const ZEvent *from, const ZEvent *to)
{
  assert(from && to && "Null pointer event");
  unsigned li = graph.lineID(from);
  unsigned lj = graph.lineID(to);
  assert(li != lj);
  unsigned li_evx = from->event_id();
  unsigned lj_evx = to->event_id();
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
       lj_after_evx < graph.lines[lj].size();
       ++lj_after_evx) {
    if (_pred[lj][li][lj_after_evx] >= (int) li_evx)
      break; // since for bigger indices also >= li_evx
    _pred[lj][li][lj_after_evx] = (int) li_evx;
  }
}


void ZPartialOrder::addLine(const ZEvent * ev)
{
  // This should be called right after line is added to graph
  assert(ev && "Null pointer event");
  assert(ev->event_id() == 0);
  assert(ev->thread_id() < 100);
  assert(graph.size() == _succ.size() + 1 &&
         graph.size() == _pred.size() + 1);
  assert(graph.lineID(ev) == _succ.size());

  // Clocks from original lines to new one
  for (unsigned li=0; li<_succ.size(); li++) {
    assert(_succ.size() == _succ[li].size());
    assert(_pred.size() == _pred[li].size());
    _succ[li].push_back(std::vector<int>(graph.lines[li].size(), INT_MAX));
    _pred[li].push_back(std::vector<int>(graph.lines[li].size(), -1));
    assert(graph.size() == _succ[li].size());
    assert(graph.size() == _pred[li].size());
  }

  // Clocks for new line
  unsigned lnew = _succ.size();
  _succ.push_back(std::vector<std::vector<int>>());
  _succ[lnew].reserve(graph.size());
  _pred.push_back(std::vector<std::vector<int>>());
  _pred[lnew].reserve(graph.size());
  assert(graph.size() == _succ.size() && graph.size() == _pred.size());
  assert(graph.lines.at(lnew).size() == 0);
  for (unsigned lold=0; lold<graph.size(); lold++) {
    _succ[lnew].push_back(std::vector<int>());
    _pred[lnew].push_back(std::vector<int>());
  }
}


void ZPartialOrder::addEvent(const ZEvent * ev)
{
  assert(ev && "Null pointer event");
  assert(graph.hasEvent(ev));
  unsigned lID = graph.lineID(ev);
  unsigned evID = ev->event_id();

  // This should be called right after ev is added to graph
  assert(graph.size() == _succ.size() && graph.size() == _pred.size());
  assert(evID == graph.lines[lID].size() - 1);
  for (unsigned li=0; li<graph.size(); li++) {
    assert(evID == _succ[lID][li].size());
    assert(evID == _pred[lID][li].size());
    // New _succ slot should have INT_MAX,
    _succ[lID][li].push_back(INT_MAX);
    // New _pred slot should have what the last slot says (if there is any; else -1)
    int newpred = (evID > 0) ? _pred[lID][li][evID-1] : -1;
    _pred[lID][li].push_back(newpred);
    assert(graph.lines[lID].size() == _succ[lID][li].size());
    assert(graph.lines[lID].size() == _pred[lID][li].size());
  }
}


void ZPartialOrder::shrink()
{
  assert(_succ.size() == _pred.size());
  _succ.shrink_to_fit();
  _pred.shrink_to_fit();
  for (unsigned i=0; i<_succ.size(); ++i) {
    assert(_succ[i].size() == _pred[i].size());
    _succ[i].shrink_to_fit();
    _pred[i].shrink_to_fit();
    for (unsigned j=0; j<_succ[i].size(); ++j) {
      assert(_succ[i][j].size() == _pred[i][j].size());
      _succ[i][j].shrink_to_fit();
      _pred[i][j].shrink_to_fit();
    }
  }
}


std::string ZPartialOrder::to_string() const
{
  std::stringstream res;

  res << "\ndigraph {\n";

  const auto& th_aux = graph.threads_auxes;

  // NODES
  for (unsigned tid = 0; tid < th_aux.size(); ++tid) {
    res << "subgraph cluster_" << tid << "{\n";
    res << "style=\"bold,rounded\" label = \"Th" << tid
        << " " << graph(tid, -1)[0]->cpid() << "\"\n";
    for (const auto& aux : th_aux[tid]) {
      unsigned line =  graph.thread_aux_to_line_id.at(std::pair<unsigned,int>(tid, aux));
      res << "subgraph cluster_" << 1001+tid*100+aux << "{\n";
      res << "style=\"invis\" label = \"" << ((aux==-1)?"real":"aux") << "\"\n";
      for (unsigned evid = 0; evid < graph(tid, aux).size(); ++evid) {
        const ZEvent *ev = graph(tid, aux)[evid];
        res << "NODE" << line * 100000 + evid
            << " [shape=\"rectangle\", label=\"" << ev->to_string(false) << "\"]\n";
      }

      res << "}\n";
    }

    res << "}\n";
  }
  res << "\n";

  // THREAD ORDER
  for (unsigned tid = 0; tid < th_aux.size(); ++tid) {
    for (const auto& aux : th_aux[tid]) {
      unsigned line =  graph.thread_aux_to_line_id.at(std::pair<unsigned,int>(tid, aux));
      for (int evid = 0; evid < (int) graph(tid, aux).size() - 1; ++evid) {
        res << "NODE" << line * 100000 + evid
            << " -> NODE" << line * 100000 + evid + 1 << "[style=bold]\n";
      }
    }
  }

  // REST ORDER
  for (unsigned tidI = 0; tidI < th_aux.size(); ++tidI) {
  for (const auto& auxI : th_aux.at(tidI)) {
    unsigned lI = graph.thread_aux_to_line_id.at(std::pair<unsigned,int>(tidI, auxI));
    unsigned realI = (auxI == -1) ? lI
    : graph.thread_aux_to_line_id.at(std::pair<unsigned,int>(tidI, -1));

    for (unsigned tidJ = 0; tidJ < th_aux.size(); ++tidJ) {
    for (const auto& auxJ : th_aux.at(tidJ)) {
      unsigned lJ = graph.thread_aux_to_line_id.at(std::pair<unsigned,int>(tidJ, auxJ));
      unsigned realJ = (auxJ == -1) ? lJ
      : graph.thread_aux_to_line_id.at(std::pair<unsigned,int>(tidJ, -1));

      if (lI != lJ) {
        int curval = INT_MAX;
        for (int liev = _succ.at(lI).at(lJ).size() - 1; liev >= 0; --liev) {
          if (auxI != -1 && realI != lJ) {
            // Do not add edges that transitively follow from
            // lI(aux) -> realI -> lJ
            auto fromev = graph.getEvent(tidI, auxI, liev);
            auto succ_idx = succ(fromev, tidI, -1);
            if (succ_idx.second < (int) _succ.at(realI).at(lJ).size() &&
                _succ.at(realI).at(lJ).at(succ_idx.second) < curval)
              curval = _succ.at(realI).at(lJ).at(succ_idx.second);
          }
          if (auxJ != -1 && lI != realJ) {
            // Do not add edges that transitively follow from
            // lI -> realJ -> lJ(aux)
            auto fromev = graph.getEvent(tidI, auxI, liev);
            auto succ_idx = succ(fromev, tidJ, -1);
            if (succ_idx.second < (int) _succ.at(realJ).at(lJ).size() &&
                _succ.at(realJ).at(lJ).at(succ_idx.second) < curval)
              curval = _succ.at(realJ).at(lJ).at(succ_idx.second);
          }
          if (_succ.at(lI).at(lJ).at(liev) < curval) {
            curval = _succ.at(lI).at(lJ).at(liev);
              res << "NODE" << lI * 100000 + liev
                  << " -> NODE" << lJ * 100000 + curval << "\n";
          }
        }
      }
    }
    }
  }
  }

  res << "}\n\n";

  return res.str();
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZPartialOrder& po)
{
  out << po.to_string();
  return out;
}


void ZPartialOrder::dump() const
{
  llvm::errs() << *this << "\n";
}
