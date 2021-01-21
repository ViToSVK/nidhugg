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
#include "ZExplorer.h"


ZGraph ZPartialOrder::graphDummy = ZGraph();


// Empty
ZPartialOrder::ZPartialOrder()
  : graph(ZPartialOrder::graphDummy),
    _succ(),
    _pred(),
    _closure_safe_until(),
    _line_sizes(),
    _threads_spanned()
{
  assert(empty());
}


// Initial
ZPartialOrder::ZPartialOrder(const ZGraph& graph)
  : graph(graph),
    _succ(),
    _pred(),
    _closure_safe_until(),
    _line_sizes(),
    _threads_spanned()
{
  assert(empty());
}


// When extending
ZPartialOrder::ZPartialOrder(ZPartialOrder&& oth, const ZGraph& graph)
  : graph(graph),
    _succ(std::move(oth._succ)),
    _pred(std::move(oth._pred)),
    _closure_safe_until(std::move(oth._closure_safe_until)),
    _line_sizes(std::move(oth._line_sizes)),
    _threads_spanned(std::move(oth._threads_spanned))
{
  // We extend from a closed partial order, hence everything is safe
  assert(graph.size() == graph._cpid_to_line.size());
  assert(graph.size() == graph._line_to_cpid.size());
  for (unsigned i = 0; i < _closure_safe_until.size(); ++i) {
    _closure_safe_until[i] = line_size(i);
  }
}


const std::vector<CPid>& ZPartialOrder::threads_spanned() const
{
  return _threads_spanned;
}


std::map<CPid, int> ZPartialOrder::thread_sizes_minus_one() const
{
  std::map<CPid, int> res;
  for (const CPid& cpid : threads_spanned()) {
    assert(spans_thread(cpid));
    assert(!res.count(cpid));
    assert(thread_size(cpid) > 0);
    res.emplace(cpid, thread_size(cpid) - 1);
  }
  return res;
}


int ZPartialOrder::line_size(unsigned line_id) const
{
  assert(line_id < _succ.size() && line_id < _line_sizes.size());
  assert(line_id < _threads_spanned.size());
  int size = _line_sizes[line_id];
  assert(size >= 0);
  return size;
}


bool ZPartialOrder::spans_thread(const CPid& cpid) const
{
  assert(graph.has_thread(cpid));
  unsigned line_id = graph.line_id(cpid);
  assert(_succ.size() == _pred.size());
  assert(_succ.size() == _threads_spanned.size());
  bool res = (line_id < _threads_spanned.size());
  assert(!res || _threads_spanned[line_id] == cpid);
  return res;
}


int ZPartialOrder::thread_size(const CPid& cpid) const
{
  assert(spans_thread(cpid));
  unsigned line_id = graph.line_id(cpid);
  return line_size(line_id);
}


bool ZPartialOrder::spans_event(const ZEvent *ev) const
{
  assert(ev);
  assert(graph.has_event(ev));
  if (ev == graph.initial())
    return true;
  if (!spans_thread(ev->cpid()))
    return false;
  return (thread_size(ev->cpid()) > ev->event_id());
}


bool ZPartialOrder::is_closure_safe(const ZEvent *read) const
{
  assert(is_read(read));
  assert(graph.has_event(read));
  assert(spans_event(read));
  assert(_closure_safe_until.size() > graph.line_id(read->cpid()));
  return (read->event_id() < _closure_safe_until[graph.line_id(read->cpid())]);
}


std::pair<const ZEvent *, int> ZPartialOrder::succ(const ZEvent *from, const CPid& to_cpid) const
{
  assert(from && "Null pointer event");
  assert(from != graph.initial() && "Asking succ for the initial node");
  assert(graph.has_event(from));
  assert(spans_event(from));
  assert(graph.has_thread(to_cpid) && "Such thread does not exist");
  assert(spans_thread(to_cpid) && "Not spanning such thread");
  assert(from->cpid() != to_cpid && "Same thread");

  unsigned from_line = graph.line_id(from);
  unsigned to_line = graph.line_id(to_cpid);
  assert(from_line != to_line);
  assert(from_line < _succ.size() && to_line < _succ.size());
  assert(to_line < _succ.at(from_line).size());
  assert(from->event_id() < _succ.at(from_line).at(to_line).size());
  int succ_evid = _succ[from_line][to_line][from->event_id()];
  assert(succ_evid >= 0);

  if (succ_evid >= line_size(to_line)) {
    assert(succ_evid == INT_MAX);
    return {nullptr, succ_evid};
  }
  assert(succ_evid < line_size(to_line));
  assert(to_line < graph._lines.size());
  assert(succ_evid < graph._lines[to_line].size());
  return {graph._lines[to_line][succ_evid], succ_evid};
}


std::pair<const ZEvent *, int> ZPartialOrder::pred(const ZEvent *to, const CPid& from_cpid) const
{
  assert(to && "Null pointer event");
  assert(to != graph.initial() && "Asking pred for the initial node");
  assert(graph.has_event(to));
  assert(spans_event(to));
  assert(graph.has_thread(from_cpid) && "Such thread does not exist");
  assert(spans_thread(from_cpid) && "Not spanning such thread");
  assert(to->cpid() != from_cpid && "Same thread");

  unsigned from_line = graph.line_id(from_cpid);
  unsigned to_line = graph.line_id(to);
  assert(from_line != to_line);
  assert(from_line < _pred.size() && to_line < _pred.size());
  assert(from_line < _pred.at(to_line).size());
  assert(to->event_id() < _pred.at(to_line).at(from_line).size());
  int pred_evid = _pred[to_line][from_line][to->event_id()];
  assert(pred_evid >= -1);

  if (pred_evid < 0) {
    assert(pred_evid == -1);
    return {nullptr, pred_evid};
  }
  assert(pred_evid < line_size(from_line));
  assert(from_line < graph._lines.size());
  assert(pred_evid < graph._lines[from_line].size());
  return {graph._lines[from_line][pred_evid], pred_evid};
}


bool ZPartialOrder::has_edge(const ZEvent *from, const ZEvent *to) const
{
  assert(from && to && "Null pointer event");
  assert(from != to && *from != *to);
  assert(graph.has_event(from) && graph.has_event(to));
  assert(spans_event(from) && spans_event(to));
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
    assert(spans_event(read_of_to));
    assert(same_ml(read_of_to, to));
    assert(!are_ordered(from, read_of_to));
    add_edge(from, read_of_to);
    return;
  }

  if (from->is_read_of_rmw() ||
      (from->is_read_of_cas() && from->value() == from->cas_compare_val())) {
    // not just 'from', but also its write-part of
    // the atomic event needs to happen before 'to'
    assert(graph.has_event(from->cpid(), from->event_id() + 1));
    // the write-part is always present in the graph, but maybe not in po
    const ZEvent * write_of_from = graph.event(from->cpid(), from->event_id() + 1);
    if (spans_event(write_of_from)) {
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

  assert(_succ.size() == _pred.size());
  unsigned li = graph.line_id(from);
  unsigned lj = graph.line_id(to);
  assert(li != lj);
  assert(li < _succ.size() && lj < _succ.size());

  std::set<const ZEvent *> before_from, after_to;
  for (unsigned lk = 0; lk<_succ.size(); ++lk) {
    if (lk != li && lk != lj) {
      auto maxbefore = pred(from, graph.line_id_to_cpid(lk));
      assert(maxbefore.second < line_size(lk));
      if (maxbefore.first)
        before_from.emplace(maxbefore.first);
      else assert(maxbefore.second == -1);

      auto minafter = succ(to, graph.line_id_to_cpid(lk));
      assert(minafter.second >= 0);
      if (minafter.first)
        after_to.emplace(minafter.first);
      else
        assert(minafter.second >= line_size(lk));
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
  assert(spans_event(from) && spans_event(to));
  unsigned li = graph.line_id(from);
  unsigned lj = graph.line_id(to);
  assert(li != lj);
  unsigned li_evx = from->event_id();
  unsigned lj_evx = to->event_id();
  assert(_succ.size() == _pred.size());
  assert(li < _succ.size() && lj < _succ.size());
  assert(_succ[li].size() == _pred[li].size());
  assert(_succ[lj].size() == _pred[lj].size());
  assert(lj < _succ[li].size() && li < _succ[lj].size());
  assert(_succ[li][lj].size() == _pred[li][lj].size());
  assert(_succ[lj][li].size() == _pred[lj][li].size());
  assert(li_evx < _succ[li][lj].size() && lj_evx < _succ[lj][li].size());
  assert( _succ[li][lj][li_evx] > (int) lj_evx && // ! li[li_evx] HB lj[lj_evx]
          _succ[lj][li][lj_evx] > (int) li_evx && // ! lj[lj_evx] HB li[li_evx]
          "Tried to add an edge between ordered nodes");
  assert( _pred[lj][li][lj_evx] < (int) li_evx && // ! li[li_evx] HB lj[lj_evx]
          _pred[li][lj][li_evx] < (int) lj_evx && // ! lj[lj_evx] HB li[li_evx]
          "Inconsistent succ/pred vector clocks");

  // Maintenance of closure_safe_until
  assert(graph.size() == graph._cpid_to_line.size());
  assert(graph.size() == graph._line_to_cpid.size());
  assert(_succ.size() == _closure_safe_until.size());
  assert(_pred.size() == _closure_safe_until.size());
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
      assert(lj < _pred.size() && lk < _pred[lj].size() &&
             lj_evx < _pred[lj][lk].size());
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
       lj_after_evx < line_size(lj);
       ++lj_after_evx) {
    assert(lj_after_evx < _pred[lj][li].size());
    if (_pred[lj][li][lj_after_evx] >= (int) li_evx)
      break; // since for bigger indices also >= li_evx
    _pred[lj][li][lj_after_evx] = (int) li_evx;
  }
}


void ZPartialOrder::add_line(const CPid& cpid)
{
  assert(graph.has_thread(cpid));
  assert(!spans_thread(cpid));
  assert(_succ.size() == _pred.size());
  _line_sizes.push_back(0);
  _threads_spanned.push_back(cpid);
  assert(_succ.size() + 1 == _line_sizes.size());
  assert(_threads_spanned.size() == _line_sizes.size());

  // Clocks from original lines to new one
  for (unsigned li=0; li<_succ.size(); li++) {
    assert(_succ.size() == _succ[li].size());
    assert(_pred.size() == _pred[li].size());
    _succ[li].push_back(std::vector<int>(line_size(li), INT_MAX));
    _pred[li].push_back(std::vector<int>(line_size(li), -1));
    assert(_succ.size() + 1 == _succ[li].size());
    assert(_pred.size() + 1 == _pred[li].size());
  }

  // Clocks for new line
  unsigned lnew = _succ.size();
  _succ.push_back(std::vector<std::vector<int>>());
  _succ[lnew].reserve(_succ.size());
  _pred.push_back(std::vector<std::vector<int>>());
  _pred[lnew].reserve(_pred.size());
  assert(_succ.size() == _pred.size() && graph.size() >= _succ.size());
  assert(line_size(lnew) == 0);
  for (unsigned lold=0; lold<_succ.size(); lold++) {
    _succ[lnew].push_back(std::vector<int>());
    _pred[lnew].push_back(std::vector<int>());
  }
  assert(_succ.size() == _line_sizes.size());

  // Nothing is closure-safe
  for (unsigned i=0; i<_closure_safe_until.size(); ++i)
    _closure_safe_until[i] = 0;
  // Closure-safe into for the new thread
  assert(_succ.size() == _closure_safe_until.size() + 1);
  _closure_safe_until.push_back(0);
  assert(_succ.size() == _closure_safe_until.size());
  assert(spans_thread(cpid));
}


void ZPartialOrder::add_line(const ZEvent * ev)
{
  assert(ev && "Null pointer event");
  assert(ev->event_id() == 0);
  assert(graph.line_id(ev) < graph._lines.size());
  assert(!graph._lines.at(graph.line_id(ev)).empty());
  assert(ev == graph._lines.at(graph.line_id(ev)).at(0));
  assert(graph.size() > _succ.size() &&
         graph.size() > _pred.size());
  assert(graph.line_id(ev) >= _succ.size());
  add_line(ev->cpid());
}


void ZPartialOrder::inherit_lines(const ZPartialOrder& oth)
{
  assert(empty());
  for (const CPid& cpid : oth.threads_spanned()) {
    add_line(cpid);
  }
}


void ZPartialOrder::add_event(const ZEvent * ev)
{
  assert(ev && "Null pointer event");
  assert(graph.has_event(ev));
  assert(!spans_event(ev));
  unsigned lID = graph.line_id(ev);
  unsigned evID = ev->event_id();

  assert(_succ.size() == _pred.size());
  assert(evID < graph._lines[lID].size());
  assert(evID == line_size(lID));
  assert(lID < _line_sizes.size());
  ++_line_sizes[lID];
  assert(evID == line_size(lID) - 1);
  for (unsigned li=0; li<_succ.size(); li++) {
    if (li == lID) {
      assert(_succ[lID][li].empty());
      assert(_pred[lID][li].empty());
      continue;
    }
    assert(evID == _succ[lID][li].size());
    assert(evID == _pred[lID][li].size());
    // New _succ slot should have INT_MAX,
    _succ[lID][li].push_back(INT_MAX);
    // New _pred slot should have what the last slot says (if there is any; else -1)
    int newpred = (evID > 0) ? _pred[lID][li][evID-1] : -1;
    _pred[lID][li].push_back(newpred);
    assert(line_size(lID) == _succ[lID][li].size());
    assert(line_size(lID) == _pred[lID][li].size());
  }
  assert(evID == line_size(lID) - 1);

  // Maintenance of closure_safe_until
  assert(graph.size() == graph._cpid_to_line.size());
  assert(graph.size() == graph._line_to_cpid.size());
  assert(_succ.size() == _closure_safe_until.size());
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
  assert(spans_event(ev));
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
  _line_sizes.shrink_to_fit();
  _threads_spanned.shrink_to_fit();
}


void ZPartialOrder::extend
(const ZEvent *read_lock, const ZAnnotation& mutated_annotation,
 const ZExplorer& explorer, const ZPartialOrder& po_full)
{
  assert((is_read(read_lock) && mutated_annotation.defines(read_lock)) ||
         (is_lock(read_lock) && mutated_annotation.is_last_lock(read_lock)));
  assert(spans_event(read_lock));

  std::list<CPid> todo;
  todo.push_back(read_lock->cpid());

  std::map<CPid, const ZEvent *> spawns;

  while (!todo.empty()) {
    CPid curr_cpid(todo.front());
    todo.pop_front();
    assert(graph.has_thread(curr_cpid));
    assert(curr_cpid == read_lock->cpid() || !spans_thread(curr_cpid) ||
           (is_join(graph.event(curr_cpid, thread_size(curr_cpid)))));

    if (!spans_thread(curr_cpid))
      add_line(curr_cpid);
    assert(spans_thread(curr_cpid));

    int full_size = graph(curr_cpid).size();
    for (int ev_id = thread_size(curr_cpid); ev_id < full_size; ++ev_id) {
      const ZEvent *ev = graph.event(curr_cpid, ev_id);
      assert(!spans_event(ev));
      assert(ev->event_id() == ev_id);

      // Handle join before adding the event
      if (is_join(ev)) {
        assert(ev->event_id() > 0);
        assert(graph.has_thread(ev->childs_cpid()));
        if (!spans_thread(ev->childs_cpid()) ||
            thread_size(ev->childs_cpid()) < graph(ev->childs_cpid()).size()) {
          // Join waiting for unfinished thread, do not add it and break
          break;
        }
        assert(spans_thread(ev->childs_cpid()) &&
               thread_size(ev->childs_cpid()) == graph(ev->childs_cpid()).size());
      }

      add_event(ev);
      assert(spans_event(ev));

      // Spawn edge
      if (ev->event_id() == 0) {
        assert(spawns.count(ev->cpid()));
        const ZEvent *spwn = spawns[ev->cpid()];
        assert(spans_event(spwn));
        assert(!are_ordered(spwn, ev));
        add_edge(spwn, ev);
      }

      // Join edge
      if (is_join(ev)) {
        assert(spans_thread(ev->childs_cpid()) &&
               thread_size(ev->childs_cpid()) == graph(ev->childs_cpid()).size());
        int lastev_idx = graph(ev->childs_cpid()).size() - 1;
        assert(lastev_idx >= 0);
        const ZEvent *wthr = graph.event(ev->childs_cpid(), lastev_idx);
        assert(spans_event(wthr));
        assert(!has_edge(ev, wthr));
        if (!has_edge(wthr, ev)) {
          add_edge(wthr, ev);
        }
      }

      // Handle specific event types
      if (is_write(ev) || is_lock(ev)) {
        if (explorer.parents.count(ev->ml())) {
          // Process backtrack points
          explorer.process_backtrack_points(po_full, ev);
        }
        if (explorer.waitfor_negallowed.count(ev->ml())) {
          // Process backtrack points
          explorer.process_backtrack_points_negallowed(po_full, ev);
        }
      }
      if (is_read(ev) || is_lock(ev)) {
        // We reached new unannotated event, break extending in curr_cpid
        break;
      }
      if (is_spawn(ev)) {
        assert(!spawns.count(ev->childs_cpid()));
        spawns.emplace(ev->childs_cpid(), ev);
        assert(graph.has_thread(ev->childs_cpid()));
        assert(!spans_thread(ev->childs_cpid()));
        todo.push_back(ev->childs_cpid());
      }
    }
    // We finished extending in curr_cpid
    if (thread_size(curr_cpid) < full_size)
      continue;
    assert(thread_size(curr_cpid) == full_size);
    // We extended until the complete end of curr_cpid
    // Locate threads waiting to join curr_cpid
    for (const CPid& jn_cpid : threads_spanned()) {
      int jn_line_size = thread_size(jn_cpid);
      assert(jn_line_size > 0);
      if (jn_line_size >= graph(jn_cpid).size()) {
        // jn_cpid is already fully extended
        assert(jn_line_size == graph(jn_cpid).size());
        continue;
      }
      const ZEvent * jn_ev = graph.event(jn_cpid, jn_line_size - 1);
      assert(spans_event(jn_ev));
      if ((is_read(jn_ev) && !mutated_annotation.defines(jn_ev)) ||
          (is_lock(jn_ev) && !mutated_annotation.is_last_lock(jn_ev))) {
        // Do not extend jn_cpid if our last event of it is unannotated
        continue;
      }
      jn_ev = graph.event(jn_cpid, jn_line_size);
      assert(!spans_event(jn_ev));
      if (is_join(jn_ev) && jn_ev->childs_cpid() == curr_cpid) {
        // Located: next event of jn_cpid wants to join curr_cpid
        todo.push_back(jn_cpid);
      }
    }
  }
}


void ZPartialOrder::process_remaining_events_for_backtrack_points
(const ZExplorer& explorer, const ZPartialOrder& po_full) const
{
  for (const auto& cpid_lineid : graph.threads()) {
    const CPid& cpid = cpid_lineid.first;
    int start = spans_thread(cpid) ? thread_size(cpid) : 0;
    for (int evidx = start; evidx < graph(cpid).size(); ++evidx) {
      assert(graph.has_event(cpid, evidx));
      const ZEvent * ev = graph.event(cpid, evidx);
      if (is_write(ev) || is_lock(ev)) {
        if (explorer.parents.count(ev->ml())) {
          // Process backtrack points
          explorer.process_backtrack_points(po_full, ev);
        }
        if (explorer.waitfor_negallowed.count(ev->ml())) {
          // Process backtrack points
          explorer.process_backtrack_points_negallowed(po_full, ev);
        }
      }
    }
  }
}


bool ZPartialOrder::something_to_annotate(const ZAnnotation& annotation) const
{
  for (const auto& cpid_lastread : graph.cache().last_read) {
    const CPid& cpid = cpid_lastread.first;
    const ZEvent * lastread = cpid_lastread.second;
    assert(is_read(lastread));
    assert(cpid == lastread->cpid());
    assert(graph.has_event(lastread));
    if (!annotation.defines(lastread)) {
      // Unannotated read
      return true;
    }
    assert(spans_event(lastread));
  }
  for (const auto& cpid_lastlock : graph.cache().last_lock) {
    const CPid& cpid = cpid_lastlock.first;
    const ZEvent * lastlock = cpid_lastlock.second;
    assert(is_lock(lastlock));
    assert(cpid == lastlock->cpid());
    assert(graph.has_event(lastlock));
    if (!spans_event(lastlock)) {
      // Lock outside of po_part
      return true;
    }
    if (lastlock->event_id() == thread_size(cpid) - 1 &&
        !annotation.is_last_lock(lastlock)) {
      // Lock last in its thread in po_part, and not annotated
      return true;
    }
  }
  return false;
}


std::string ZPartialOrder::to_string() const
{
  std::stringstream res;
  assert(_succ.size() == _pred.size());

  res << "\ndigraph {\n";

  // NODES
  for (unsigned tid = 0; tid < _succ.size(); ++tid) {
    res << "subgraph cluster_" << tid << "{\n";
    res << "style=\"bold,rounded\" label = \"Th" << tid
        << " " << graph.line_id_to_cpid(tid) << "\"\n";
    for (unsigned evid = 0; evid < line_size(tid); ++evid) {
      const ZEvent *ev = graph._lines[tid][evid];
      res << "NODE" << tid * 100000 + evid
          << " [shape=\"rectangle\", label=\"" << ev->to_string(false) << "\"]\n";
    }
    res << "}\n";
  }
  res << "\n";

  // THREAD ORDER
  for (unsigned tid = 0; tid < _succ.size(); ++tid) {
    // evid is int to prevent unsigned underflow of the for-upper-bound
    for (int evid = 0; evid < line_size(tid) - 1; ++evid) {
      res << "NODE" << tid * 100000 + evid
          << " -> NODE" << tid * 100000 + evid + 1 << "[style=bold]\n";
    }
  }

  // REST ORDER
  for (unsigned ti = 0; ti < _succ.size(); ++ti) {
    for (unsigned tj = 0; tj < _succ.size(); ++tj) {
      if (ti != tj) {
        int curval = INT_MAX;
        for (int evi = _succ.at(ti).at(tj).size() - 1; evi >= 0; --evi) {
          if (_succ.at(ti).at(tj).at(evi) < curval) {
            curval = _succ.at(ti).at(tj).at(evi);
            assert(curval < line_size(tj));
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
