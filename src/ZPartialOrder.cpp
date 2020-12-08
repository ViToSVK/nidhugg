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


// Empty
ZPartialOrder::ZPartialOrder()
  : graph(ZPartialOrder::graphDummy),
    _succ(),
    _pred(),
    _closure_safe_until()
{
  assert(empty());
}


// Initial
ZPartialOrder::ZPartialOrder(const ZGraph& graph)
  : graph(graph),
    _succ(),
    _pred(),
    _closure_safe_until()
{
  assert(empty());
}


// When extending
ZPartialOrder::ZPartialOrder(ZPartialOrder&& oth, const ZGraph& graph)
  : graph(graph),
    _succ(std::move(oth._succ)),
    _pred(std::move(oth._pred)),
    _closure_safe_until(std::move(oth._closure_safe_until))
{
  // We extend from a closed partial order, hence everything is safe
  assert(graph.size() == graph._cpid_to_line.size());
  assert(graph.size() == graph._line_to_cpid.size());
  assert(graph.size() == _closure_safe_until.size());
  for (unsigned i = 0; i < graph.size(); ++i) {
    _closure_safe_until[i] = graph._lines[i].size();
  }
}


bool ZPartialOrder::is_closure_safe(const ZEvent *read) const
{
  assert(is_read(read));
  assert(graph.has_event(read));
  assert(_closure_safe_until.size() > graph.line_id(read->cpid()));
  return (read->event_id() < _closure_safe_until[graph.line_id(read->cpid())]);
}


std::pair<const ZEvent *, int> ZPartialOrder::succ(const ZEvent *from, const CPid& to_cpid) const
{
  assert(from && "Null pointer event");
  assert(graph.has_event(from));
  assert(from != graph.initial() && "Asking succ for the initial node");
  assert(graph.has_thread(to_cpid) && "Such thread does not exist");
  assert(from->cpid() != to_cpid && "Same thread");

  unsigned from_line = graph.line_id(from);
  unsigned to_line = graph.line_id(to_cpid);
  assert(from_line != to_line);
  assert(from->event_id() < _succ.at(from_line).at(to_line).size());
  int succ_evid = _succ[from_line][to_line][from->event_id()];
  assert(succ_evid >= 0);

  assert(from_line < graph._lines.size() && to_line < graph._lines.size());
  if (succ_evid >= (int) graph._lines[to_line].size()) {
    assert(succ_evid == INT_MAX);
    return {nullptr, succ_evid};
  }
  assert(succ_evid < (int) graph._lines[to_line].size());
  return {graph._lines[to_line][succ_evid], succ_evid};
}


std::pair<const ZEvent *, int> ZPartialOrder::pred(const ZEvent *to, const CPid& from_cpid) const
{
  assert(to && "Null pointer event");
  assert(graph.has_event(to));
  assert(to != graph.initial() && "Asking pred for the initial node");
  assert(graph.has_thread(from_cpid) && "Such thread does not exist");
  assert(to->cpid() != from_cpid && "Same thread");

  unsigned from_line = graph.line_id(from_cpid);
  unsigned to_line = graph.line_id(to);
  assert(from_line != to_line);
  assert(to->event_id() < _pred.at(to_line).at(from_line).size());
  int pred_evid = _pred[to_line][from_line][to->event_id()];
  assert(pred_evid >= -1);

  assert(from_line < graph._lines.size() && to_line < graph._lines.size());
  if (pred_evid < 0) {
    assert(pred_evid == -1);
    return {nullptr, pred_evid};
  }
  assert(pred_evid < (int) graph._lines[from_line].size());
  return {graph._lines[from_line][pred_evid], pred_evid};
}


bool ZPartialOrder::has_edge(const ZEvent *from, const ZEvent *to) const
{
  assert(from && to && "Null pointer event");
  assert(from != to && *from != *to);
  assert(graph.has_event(from) && graph.has_event(to));
  if (from == graph.initial())
    return true; // init HB everything
  if (to == graph.initial())
    return false; // nothing HB init
  if (from->cpid() == to->cpid())
    return from->event_id() < to->event_id();

  unsigned from_line = graph.line_id(from);
  unsigned to_line = graph.line_id(to);
  return _succ[from_line][to_line]
              [from->event_id()] <= to->event_id();
}


bool ZPartialOrder::are_ordered(const ZEvent *ev1, const ZEvent *ev2) const
{
  return has_edge(ev1, ev2) || has_edge(ev2, ev1);
}


// Adds an edge between two unordered nodes
// This method maintains:
// 1) line-pair-wise transitivity
// 2) complete transitivity
void ZPartialOrder::add_edge(const ZEvent *from, const ZEvent *to)
{
  assert(from && to && "Null pointer event");
  assert(!are_ordered(from, to));

  if (to->is_write_of_atomic_event()) {
    // 'from' needs to happen already before
    // the read-part of the 'to' atomic event
    assert(graph.has_event(to->cpid(), to->event_id() - 1));
    const ZEvent * read_of_to = graph.event(to->cpid(), to->event_id() - 1);
    assert(read_of_to && read_of_to->is_read_of_atomic_event());
    assert(same_ml(read_of_to, to));
    assert(!are_ordered(from, read_of_to));
    add_edge(from, read_of_to);
    return;
  }

  if (from->is_read_of_atomic_event()) {
    // not just 'from', but also its write-part of
    // the atomic event needs to happen before 'to'
    if (graph.has_event(from->cpid(), from->event_id() + 1)) {
      // the write-part is present in the graph
      const ZEvent * write_of_from = graph.event(from->cpid(), from->event_id() + 1);
      assert(write_of_from && write_of_from->is_write_of_atomic_event());
      assert(same_ml(from, write_of_from));
      assert(!are_ordered(write_of_from, to));
      add_edge(write_of_from, to);
      return;
    }
  }

  // Maintenance of the complete transitivity
  // Collect nodes from different lines with edges:
  // 1) to   li[li_evx]
  // 2) from lj[lj_evx]

  assert(graph.size() == _succ.size() && graph.size() == _pred.size());
  unsigned li = graph.line_id(from);
  unsigned lj = graph.line_id(to);
  assert(li != lj);

  std::set<const ZEvent *> before_from, after_to;
  for (unsigned lk = 0; lk<graph.size(); ++lk) {
    if (lk != li && lk != lj) {
      auto maxbefore = pred(from, graph.line_id_to_cpid(lk));
      assert(maxbefore.second < (int) graph._lines[lk].size());
      if (maxbefore.first)
        before_from.emplace(maxbefore.first);
      else assert(maxbefore.second == -1);

      auto minafter = succ(to, graph.line_id_to_cpid(lk));
      assert(minafter.second >= 0);
      if (minafter.first)
        after_to.emplace(minafter.first);
      else
        assert(minafter.second >= (int) graph._lines[lk].size());
    }
  }

  add_edge_help(from, to);

  // TryAdd edges between each of *before_from* and *to*
  for (const auto& bef : before_from) {
    assert(!has_edge(to, bef) && "Cycle");
    if (!has_edge(bef, to))
      add_edge_help(bef, to);
  }

  // TryAdd edges between *from* and each of *after_to*
  for (const auto& aft : after_to) {
    assert(!has_edge(aft, from) && "Cycle");
    if (!has_edge(from, aft))
      add_edge_help(from, aft);
  }

  // TryAdd edges between each of *before_from* and each of *after_to*
  // (if they belong to different threads)
  for (auto& bef : before_from)
    for (auto& aft : after_to)
      if (graph.line_id(bef) != graph.line_id(aft)) {
        assert(!has_edge(aft, bef) && "Cycle");
        if (!has_edge(bef, aft))
          add_edge_help(bef, aft);
      }
}


// Helper method for add_edge,
// maintains line-pair-wise transitivity
void ZPartialOrder::add_edge_help(const ZEvent *from, const ZEvent *to)
{
  assert(from && to && "Null pointer event");
  unsigned li = graph.line_id(from);
  unsigned lj = graph.line_id(to);
  assert(li != lj);
  unsigned li_evx = from->event_id();
  unsigned lj_evx = to->event_id();
  assert( _succ[li][lj][li_evx] > (int) lj_evx && // ! li[li_evx] HB lj[lj_evx]
          _succ[lj][li][lj_evx] > (int) li_evx && // ! lj[lj_evx] HB li[li_evx]
          "Tried to add an edge between ordered nodes");
  assert( _pred[lj][li][lj_evx] < (int) li_evx && // ! li[li_evx] HB lj[lj_evx]
          _pred[li][lj][li_evx] < (int) lj_evx && // ! lj[lj_evx] HB li[li_evx]
          "Inconsistent succ/pred vector clocks");

  // Maintenance of closure_safe_until
  assert(graph.size() == graph._cpid_to_line.size());
  assert(graph.size() == graph._line_to_cpid.size());
  assert(graph.size() == _closure_safe_until.size());
  for (unsigned lk=0; lk<_closure_safe_until.size(); ++lk) {
    if (lk == lj) {
      // The new edge leads to this thread
      // Everything up until (not including) 'to' remains safe
      if (_closure_safe_until[lk] > to->event_id())
        _closure_safe_until[lk] = to->event_id();
    } else {
      // The new edge does not lead to this thread
      // Everything that was already happening before
      // 'to' remains safe
      int already = _pred[lj][lk][lj_evx] + 1;
      assert(already >= 0);
      if (_closure_safe_until[lk] > already)
        _closure_safe_until[lk] = already;
    }
  }

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
       lj_after_evx < graph._lines[lj].size();
       ++lj_after_evx) {
    if (_pred[lj][li][lj_after_evx] >= (int) li_evx)
      break; // since for bigger indices also >= li_evx
    _pred[lj][li][lj_after_evx] = (int) li_evx;
  }
}


void ZPartialOrder::add_line(const ZEvent * ev)
{
  // This should be called right after line is added to graph
  assert(ev && "Null pointer event");
  assert(ev->event_id() == 0);
  assert(graph.size() == _succ.size() + 1 &&
         graph.size() == _pred.size() + 1);
  assert(graph.line_id(ev) == _succ.size());

  // Clocks from original lines to new one
  for (unsigned li=0; li<_succ.size(); li++) {
    assert(_succ.size() == _succ[li].size());
    assert(_pred.size() == _pred[li].size());
    _succ[li].push_back(std::vector<int>(graph._lines[li].size(), INT_MAX));
    _pred[li].push_back(std::vector<int>(graph._lines[li].size(), -1));
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
  assert(graph._lines.at(lnew).size() == 0);
  for (unsigned lold=0; lold<graph.size(); lold++) {
    _succ[lnew].push_back(std::vector<int>());
    _pred[lnew].push_back(std::vector<int>());
  }

  // Nothing is closure-safe
  for (unsigned i=0; i<_closure_safe_until.size(); ++i)
    _closure_safe_until[i] = 0;
  // Closure-safe into for the new thread
  assert(graph.size() == _closure_safe_until.size() + 1);
  _closure_safe_until.push_back(0);
  assert(graph.size() == _closure_safe_until.size());
}


void ZPartialOrder::add_event(const ZEvent * ev)
{
  assert(ev && "Null pointer event");
  assert(graph.has_event(ev));
  unsigned lID = graph.line_id(ev);
  unsigned evID = ev->event_id();

  // This should be called right after ev is added to graph
  assert(graph.size() == _succ.size() && graph.size() == _pred.size());
  assert(evID == graph._lines[lID].size() - 1);
  for (unsigned li=0; li<graph.size(); li++) {
    assert(evID == _succ[lID][li].size());
    assert(evID == _pred[lID][li].size());
    // New _succ slot should have INT_MAX,
    _succ[lID][li].push_back(INT_MAX);
    // New _pred slot should have what the last slot says (if there is any; else -1)
    int newpred = (evID > 0) ? _pred[lID][li][evID-1] : -1;
    _pred[lID][li].push_back(newpred);
    assert(graph._lines[lID].size() == _succ[lID][li].size());
    assert(graph._lines[lID].size() == _pred[lID][li].size());
  }

  // Maintenance of closure_safe_until
  assert(graph.size() == graph._cpid_to_line.size());
  assert(graph.size() == graph._line_to_cpid.size());
  assert(graph.size() == _closure_safe_until.size());
  for (unsigned lk=0; lk<_closure_safe_until.size(); ++lk) {
    // Closure-safety doesn't change in the thread of the new event
    if (lk != lID) {
      // The new event is not in this thread
      // Everything in this thread which was already
      // happening before the new event 'ev' remains safe
      assert(_pred[lID][lk].size() - 1 == evID);
      assert(evID == 0 || _pred[lID][lk][evID] == _pred[lID][lk][evID - 1]);
      int already = _pred[lID][lk][evID] + 1;
      assert(already >= 0);
      if (_closure_safe_until[lk] > already)
        _closure_safe_until[lk] = already;
    }
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
  assert(_closure_safe_until.size() == _succ.size());
  _closure_safe_until.shrink_to_fit();
}


std::string ZPartialOrder::to_string() const
{
  std::stringstream res;

  res << "\ndigraph {\n";

  // NODES
  for (unsigned tid = 0; tid < graph.size(); ++tid) {
    res << "subgraph cluster_" << tid << "{\n";
    res << "style=\"bold,rounded\" label = \"Th" << tid
        << " " << graph.line_id_to_cpid(tid) << "\"\n";
    for (unsigned evid = 0; evid < graph._lines[tid].size(); ++evid) {
      const ZEvent *ev = graph._lines[tid][evid];
      res << "NODE" << tid * 100000 + evid
          << " [shape=\"rectangle\", label=\"" << ev->to_string(false) << "\"]\n";
    }
    res << "}\n";
  }
  res << "\n";

  // THREAD ORDER
  for (unsigned tid = 0; tid < graph.size(); ++tid) {
    // evid is int to prevent unsigned underflow of the for-upper-bound
    for (int evid = 0; evid < (int) graph._lines[tid].size() - 1; ++evid) {
      res << "NODE" << tid * 100000 + evid
          << " -> NODE" << tid * 100000 + evid + 1 << "[style=bold]\n";
    }
  }

  // REST ORDER
  for (unsigned ti = 0; ti < graph.size(); ++ti) {
    for (unsigned tj = 0; tj < graph.size(); ++tj) {
      if (ti != tj) {
        int curval = INT_MAX;
        for (int evi = _succ.at(ti).at(tj).size() - 1; evi >= 0; --evi) {
          if (_succ.at(ti).at(tj).at(evi) < curval) {
            curval = _succ.at(ti).at(tj).at(evi);
            assert(curval < graph._lines[tj].size());
              res << "NODE" << ti * 100000 + evi
                  << " -> NODE" << tj * 100000 + curval << "\n";
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
  llvm::errs() << to_string() << "\n";
}
