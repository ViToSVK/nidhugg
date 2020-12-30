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

#include "ZHelpers.h"
#include "ZGraph.h"
#include "ZBuilderSC.h"

//static const bool DEBUG = false;
//#include "ZDebug.h"


// Empty
ZGraph::ZGraph()
  : _init(ZEvent(true)),
    _lines(),
    _cpid_to_line(),
    _line_to_cpid(),
    _po(*this),
    _cache()
{
  assert(empty());
}


// Initial
ZGraph::ZGraph(const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& trace)
  : ZGraph()
{
  assert(trace);
  trace_to_po(*trace, nullptr);
}


// Extending
// Partial order that will be moved as original
// Trace and annotation that will extend this copy of the graph
ZGraph::ZGraph
(const ZGraph& oth,
 ZPartialOrder&& po,
 const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& trace,
 const ZAnnotation& annotation)
  : _init(ZEvent(true)),
    _lines(oth._lines),
    _cpid_to_line(oth._cpid_to_line),
    _line_to_cpid(oth._line_to_cpid),
    _po(std::move(po), *this),
    _cache()
{
  assert(trace);
  trace_to_po(*trace, &annotation);
}


/* *************************** */
/* LINES                       */
/* *************************** */

const LineT& ZGraph::operator()(const CPid& cpid) const
{
  assert(has_thread(cpid));
  auto it = _cpid_to_line.find(cpid);
  assert(it != _cpid_to_line.end());
  assert(it->second < _lines.size());
  return _lines[it->second];
}


const ZEvent * ZGraph::event(const CPid& cpid, int event_id) const
{
  assert(has_thread(cpid));
  if (event_id < 0) {
    assert(event_id == -1);
    return initial();
  }
  const LineT& line = this->operator()(cpid);
  assert(event_id < line.size());
  return line[event_id];
}


const ZEvent * ZGraph::event(const ZEventID& id) const
{
  return event(id.cpid(), id.event_id());
}


const ZEvent * ZGraph::unlock_of_this_lock(const ZEventID& id) const
{
  int cur_ev = id.event_id();
  const ZEvent * lock = event(id);
  assert(is_lock(lock));
  while (cur_ev + 1 < this->operator()(id.cpid()).size()) {
    cur_ev++;
    const ZEvent * res = event(id.cpid(), cur_ev);
    if (is_unlock(res) && same_ml(lock, res))
      return res;
  }
  return nullptr;
}


void ZGraph::add_line(const ZEvent *ev)
{
  assert(ev && !is_initial(ev));
  assert(!has_event(ev));
  assert(!has_thread(ev->cpid()));

  _cpid_to_line.emplace(ev->cpid(), _lines.size());
  _line_to_cpid.push_back(ev->cpid());
  assert(_cpid_to_line.size() == _line_to_cpid.size());
  _lines.push_back(LineT());
  _lines.back().reserve(16);
}


void ZGraph::add_event(const ZEvent *ev)
{
  assert(ev && !is_initial(ev));
  assert(!has_event(ev));
  assert(has_thread(ev->cpid()));

  auto it = _cpid_to_line.find(ev->cpid());
  assert(it != _cpid_to_line.end());
  LineT& line = _lines[it->second];

  int event_id = line.size();
  assert(event_id == ev->event_id());
  line.push_back(ev);
  assert(has_event(ev));
}


void ZGraph::replace_event(const ZEvent *old_ev, const ZEvent *new_ev)
{
  assert(old_ev && new_ev);
  assert(*old_ev == *new_ev);
  assert(old_ev->kind == new_ev->kind &&
         old_ev->id() == new_ev->id() &&
         old_ev->ml() == new_ev->ml());
  if (old_ev == new_ev) {
    // Latest change: we can have same pointers across runs now.
    return;
  }
  // assert(old_ev != new_ev && "different pointers but same ZEventID is expected");
  assert(has_event(old_ev) && !has_event(new_ev));

  unsigned lid = line_id(old_ev->cpid());
  assert((*this)(old_ev->cpid()).at(old_ev->event_id()) == old_ev);
  assert(_lines.at(lid).at(old_ev->event_id()) == old_ev);
  _lines[lid][old_ev->event_id()] = new_ev;
  assert((*this)(old_ev->cpid()).at(old_ev->event_id()) == new_ev);
  assert(_lines.at(lid).at(old_ev->event_id()) == new_ev);
  assert(!has_event(old_ev) && has_event(new_ev));
}


bool ZGraph::has_event(const CPid& cpid, int event_id) const
{
  if (!has_thread(cpid))
    return false;
  const LineT& line = this->operator()(cpid);
  assert(event_id >= 0);
  return (event_id < line.size());
}


bool ZGraph::has_event(const ZEvent *ev) const
{
  assert(ev);
  if (ev == initial())
    return true;
  if (!has_thread(ev->cpid()))
    return false;
  const LineT& line = this->operator()(ev->cpid());
  assert(ev->event_id() >= 0);
  return (ev->event_id() < line.size() && line[ev->event_id()] == ev);
}


std::map<CPid, int> ZGraph::thread_sizes_minus_one() const
{
  std::map<CPid, int> res;
  for (unsigned thr = 0; thr < _lines.size(); ++thr) {
    assert(!_lines[thr].empty());
    CPid cpid = line_id_to_cpid(thr);
    assert(!res.count(cpid));
    res.emplace(cpid, _lines[thr].size() - 1);
  }
  return res;
}


void ZGraph::shrink()
{
  _lines.shrink_to_fit();
  for (unsigned i=0; i<_lines.size(); ++i)
    _lines[i].shrink_to_fit();
  _line_to_cpid.shrink_to_fit();
  for (const auto& ml_x : _cache.writes)
    for (const auto& cpid_w : ml_x.second)
      _cache.writes[ml_x.first][cpid_w.first].shrink_to_fit();
}


unsigned ZGraph::line_id(const CPid& cpid) const
{
  assert(has_thread(cpid));
  auto it = _cpid_to_line.find(cpid);
  assert(it != _cpid_to_line.end());
  return it->second;
}


unsigned ZGraph::line_id(const ZEvent *ev) const
{
  assert(ev);
  assert(!is_initial(ev) && "Called for initial event");
  return line_id(ev->cpid());
}


const CPid& ZGraph::line_id_to_cpid(unsigned line_id) const
{
  assert(line_id < _lines.size());
  assert(line_id < _line_to_cpid.size());
  return _line_to_cpid[line_id];
}


std::set<CPid> ZGraph::threads() const
{
  std::set<CPid> res;
  for (const auto& cpid_lineid : _cpid_to_line) {
    assert(!res.count(cpid_lineid.first));
    res.insert(cpid_lineid.first);
  }
  return res;
}


bool ZGraph::has_thread(const CPid& cpid) const
{
  return _cpid_to_line.count(cpid);
}


/* *************************** */
/* GRAPH EXTENSION             */
/* *************************** */

// Extends this graph so it corresponds to 'trace'
// Check the header file for the method description
void ZGraph::trace_to_po(const std::vector<std::unique_ptr<ZEvent>>& trace,
                         const ZAnnotation *annotation_ptr)
{
  assert(!trace.empty());
  assert(trace.size() >= events_size());

  std::unordered_map<CPid, unsigned> cur_evidx;
  cur_evidx.reserve(8);
  for (const CPid& cpid : threads())
    cur_evidx.emplace(cpid, 0);

  std::unordered_set<CPid> has_unannotated_read_or_lock;
  has_unannotated_read_or_lock.reserve(8);

  std::unordered_set<SymAddrSize> found_last_lock_for_location;
  found_last_lock_for_location.reserve(8);

  std::unordered_set<CPid> forbidden_threads;
  forbidden_threads.reserve(8);

  std::unordered_set<CPid> allowed_threads;
  allowed_threads.reserve(8);
  allowed_threads.insert(trace[0]->cpid());

  std::unordered_set<CPid> threads_with_event_we_dont_add;
  threads_with_event_we_dont_add.reserve(8);

  std::unordered_map
    <CPid, std::unordered_map
     <SymAddrSize, const ZEvent *>> last_write;

  std::vector<const ZEvent *> spawns;
  spawns.reserve(8);

  std::vector<const ZEvent *> joins;
  joins.reserve(8);

  std::unordered_map<SymAddrSize, const ZEvent *> mutex_inits;
  mutex_inits.reserve(8);

  std::unordered_map<SymAddrSize, const ZEvent *> mutex_destroys;
  mutex_destroys.reserve(8);

  std::unordered_map
    <SymAddrSize, std::unordered_map<CPid, const ZEvent *>> mutex_first;
  mutex_first.reserve(8);

  std::unordered_map
    <SymAddrSize, std::unordered_map<CPid, const ZEvent *>> mutex_last;
  mutex_last.reserve(8);

  for (auto traceit = trace.begin(); traceit != trace.end(); ++traceit) {
    const ZEvent *ev = traceit->get();
    assert(!(ev->cpid().is_auxiliary()) && "Only real threads");

    if (forbidden_threads.count(ev->cpid())) {
      // This is a forbidden thread because it
      // happens after a so-far unannotated node
      continue;
    }
    if (!allowed_threads.count(ev->cpid())) {
      // This is a forbidden thread because it
      // is created after a so-far unannotated node
      forbidden_threads.insert(ev->cpid());
      continue;
    }
    assert(allowed_threads.count(ev->cpid()));

    // Check possible cases why we cannot add this event
    //
    // Check if we already haven't added a thread-predecessor of this event
    if (threads_with_event_we_dont_add.count(ev->cpid())) {
      assert(!has_thread(ev->cpid()) ||
             cur_evidx.at(ev->cpid()) == (*this)(ev->cpid()).size());
      continue;
    }
    // Check if thread already has an unannotated read or lock
    if (has_unannotated_read_or_lock.count(ev->cpid())) {
      // This event happens after something we need to
      // annotate first, so we don't add it into the graph
      assert(!has_thread(ev->cpid()) ||
             cur_evidx.at(ev->cpid()) == (*this)(ev->cpid()).size());
      threads_with_event_we_dont_add.insert(ev->cpid());
      continue;
    }
    // Joins of threads with an event we didn't include also
    // can not be included and neither any of their successors
    if (is_join(ev)) {
      bool dontInclude = false;
      if (!allowed_threads.count(ev->childs_cpid())) {
        dontInclude = true;
      } else {
        assert(has_thread(ev->childs_cpid()));
        if (threads_with_event_we_dont_add.count(ev->childs_cpid()))
          dontInclude = true;
      }
      if (dontInclude) {
        assert(!has_thread(ev->cpid()) ||
               cur_evidx.at(ev->cpid()) == (*this)(ev->cpid()).size());
        threads_with_event_we_dont_add.insert(ev->cpid());
        continue;
      }
    }

    // We will add the event
    // First, add line if not already present
    if (!has_thread(ev->cpid())) {
      // A new (allowed) thread appeared
      assert(allowed_threads.count(ev->cpid()));
      assert(!cur_evidx.count(ev->cpid()));
      cur_evidx.emplace(ev->cpid(), 0);

      add_line(ev);
      _po.add_line(ev);
      assert(has_thread(ev->cpid()));
    }
    assert(cur_evidx.count(ev->cpid()));
    unsigned ev_idx = cur_evidx[ev->cpid()];
    assert(ev->event_id() == (int) ev_idx);

    if (ev_idx == (*this)(ev->cpid()).size()) {
      // New event
      add_event(ev);
      _po.add_event(ev);

      // XXX: function pointer calls not handled
      if (is_spawn(ev))
        spawns.push_back(ev);
      if (is_join(ev))
        joins.push_back(ev);

    } else {
      // Already known event
      const ZEvent *oldev = event(ev->cpid(), ev->event_id());
      replace_event(oldev, ev);
    }

    ++cur_evidx[ev->cpid()];

    // Handle Spawn
    //
    if (is_spawn(ev)) {
      allowed_threads.insert(ev->childs_cpid());
      assert(size() > 0);
    }

    // Handle Write
    //
    if (is_write(ev)) {
      // Last write
      if (!last_write.count(ev->cpid()))
        last_write.emplace
          (ev->cpid(), std::unordered_map<SymAddrSize, const ZEvent *>());
      last_write[ev->cpid()][ev->ml()] = ev;
      // Cache - writes
      if (!_cache.writes.count(ev->ml()))
        _cache.writes.emplace
          (ev->ml(), std::unordered_map<CPid, std::vector<const ZEvent *>>());
      if (!_cache.writes[ev->ml()].count(ev->cpid()))
        _cache.writes[ev->ml()].emplace(ev->cpid(), std::vector<const ZEvent *>());
      _cache.writes[ev->ml()][ev->cpid()].push_back(ev);
    }

    // Handle Read
    //
    if (is_read(ev)) {
      // Check if ev is not annotated yet
      if (!annotation_ptr || !(annotation_ptr->defines(ev))) {
        // This will be the last node for the corresponding thread
        has_unannotated_read_or_lock.insert(ev->cpid());
      } else {
        assert(annotation_ptr && annotation_ptr->defines(ev));
      }
      // Cache - writes
      if (!_cache.writes.count(ev->ml()))
        _cache.writes.emplace
          (ev->ml(), std::unordered_map<CPid, std::vector<const ZEvent *>>());
      // Cache - local_write
      assert(!_cache.local_write.count(ev));
      if (last_write.count(ev->cpid()) &&
          last_write[ev->cpid()].count(ev->ml())) {
        assert(is_write(last_write.at(ev->cpid()).at(ev->ml())));
        _cache.local_write.emplace(ev, last_write[ev->cpid()][ev->ml()]);
      }
      else
        _cache.local_write.emplace(ev, nullptr);
    }

    // Handle Mutex events
    //
    if (is_mutex_init(ev)) {
      assert(!mutex_inits.count(ev->ml()));
      mutex_inits.emplace(ev->ml(), ev);
    }
    if (is_mutex_destroy(ev)) {
      assert(!mutex_destroys.count(ev->ml()));
      mutex_destroys.emplace(ev->ml(), ev);
    }
    if (is_lock(ev)) {
      // Check if ev is not annotated yet
      if (!annotation_ptr || !(annotation_ptr->location_has_some_lock(ev))) {
        // No annotated lock has held this location yet
        // This will be the last node for the corresponding thread
        has_unannotated_read_or_lock.insert(ev->cpid());
      } else {
        // Some lock has already happened on this location
        if (annotation_ptr->is_last_lock(ev)) {
          // This is the last annotated lock for this location
          assert(!found_last_lock_for_location.count(ev->ml()));
          found_last_lock_for_location.insert(ev->ml());
        } else if (found_last_lock_for_location.count(ev->ml())) {
          // This is a lock after the last annotated lock for this location
          // This will be the last node for the corresponding thread
          has_unannotated_read_or_lock.insert(ev->cpid());
        }
      }
    }
    if (is_lock(ev) || is_unlock(ev)) {
      // For each ml, for each thread, remember first access
      if (!mutex_first.count(ev->ml())) {
        mutex_first.emplace(ev->ml(), std::unordered_map<CPid, const ZEvent *>());
        mutex_first[ev->ml()].reserve(8);
      }
      if (!mutex_first[ev->ml()].count(ev->cpid()))
        mutex_first[ev->ml()].emplace(ev->cpid(), ev);
      // For each ml, for each thread, remember last access
      if (!mutex_last.count(ev->ml())) {
        mutex_last.emplace(ev->ml(), std::unordered_map<CPid, const ZEvent *>());
        mutex_last[ev->ml()].reserve(8);
      }
      mutex_last[ev->ml()][ev->cpid()] = ev;
    }

  } // end of loop for traversing trace and creating nodes

  shrink();
  _po.shrink();

  #ifndef NDEBUG
  assert(size() == cur_evidx.size());
  for (auto& thread_numevents : cur_evidx) {
    assert(thread_numevents.second == (*this)(thread_numevents.first).size()
           && "Didn't go through entire original part of the graph");
  }
  #endif

  // EDGES - spawns
  for (const ZEvent *spwn : spawns) {
    assert(has_thread(spwn->childs_cpid()));
    assert(!(*this)(spwn->childs_cpid()).empty());
    const ZEvent *nthr = event(spwn->childs_cpid(), 0);

    assert(!_po.has_edge(nthr, spwn));
    if (!_po.has_edge(spwn, nthr))
      _po.add_edge(spwn, nthr);
  }

  // EDGES - joins
  for (const ZEvent *jn : joins) {
    assert(has_thread(jn->childs_cpid()));
    assert(!(*this)(jn->childs_cpid()).empty());
    int lastev_idx = (*this)(jn->childs_cpid()).size() - 1;
    assert(lastev_idx >= 0);
    const ZEvent *wthr = event(jn->childs_cpid(), lastev_idx);\

    assert(!_po.has_edge(jn, wthr));
    if (!_po.has_edge(wthr, jn))
      _po.add_edge(wthr, jn);
  }

  // EDGES - mutex inits
  for (auto& in : mutex_inits) {
    const SymAddrSize& loc_init = in.first;
    const ZEvent *ev_init = in.second;

    if (mutex_first.count(loc_init)) {
      for (auto& cpid_ev_first : mutex_first[loc_init]) {
        const ZEvent *ev_first = cpid_ev_first.second;

        assert(!_po.has_edge(ev_first, ev_init));
        if (!_po.has_edge(ev_init, ev_first))
          _po.add_edge(ev_init, ev_first);
      }
    } else assert(!mutex_last.count(loc_init));
  }

  // EDGES - mutex destroys
  for (auto& de : mutex_destroys) {
    const SymAddrSize& loc_destroy = de.first;
    const ZEvent *ev_destroy = de.second;

    if (mutex_last.count(loc_destroy)) {
      for (auto& cpid_ev_last : mutex_last[loc_destroy]) {
        const ZEvent *ev_last = cpid_ev_last.second;

        assert(!_po.has_edge(ev_destroy, ev_last));
        if (!_po.has_edge(ev_last, ev_destroy))
          _po.add_edge(ev_last, ev_destroy);
      }
    } else assert(!mutex_first.count(loc_destroy));
  }
}


/* *************************** */
/* MAIN ALGORITHM              */
/* *************************** */

int ZGraph::get_tailw_index
(const SymAddrSize& ml, const CPid& cpid, int evX) const
{
  assert(_cache.writes.count(ml));
  assert(has_thread(cpid));

  if (evX < 0)
    return -1;
  if (!_cache.writes.at(ml).count(cpid)) // using 'at' because this method is const
    return -1;

  assert(_cache.writes.at(ml).count(cpid));
  const std::vector<const ZEvent *>& writes = _cache.writes.at(ml).at(cpid);
  if (writes.empty())
    return -1;

  unsigned ev = evX;
  int low = 0;
  if (ev < writes[low]->event_id())
    return -1;

  int high = writes.size() - 1;
  if (ev >= writes[high]->event_id())
    return high;

  assert(low < high);
  if (low + 1 == high) {
    assert(ev >= writes[low]->event_id() &&
           ev < writes[high]->event_id());
    return low;
  }

  // Low represents something that
  // can possibly be the answer
  // High represents something that
  // is above the answer
  // Do binary search
  while (true) {
    assert(low + 1 < high);
    assert(low >= 0 && high < writes.size());
    assert(ev >= writes[low]->event_id() &&
           ev < writes[high]->event_id());
    int mid = ((high - low) / 2) + low;
    assert(low < mid && mid < high);

    if (ev >= writes[mid]->event_id())
      low = mid;
    else
      high = mid;

    if (low + 1 == high) {
      assert(ev >= writes[low]->event_id() &&
             ev < writes[high]->event_id());
      return low;
    }
  }
}


const ZEvent * ZGraph::get_tailw
(const SymAddrSize& ml, const CPid& cpid, int evX) const
{
  assert(_cache.writes.count(ml));
  assert(has_thread(cpid));

  int idx = get_tailw_index(ml, cpid, evX);
  if (idx == -1)
    return nullptr;

  assert(idx < (int) _cache.writes.at(ml).at(cpid).size());
  auto res = _cache.writes.at(ml).at(cpid)[idx]; // using 'at' because this method is const
  assert(is_write(res) && res->ml() == ml &&
         res->cpid() == cpid &&
         evX >= 0 && res->event_id() <= evX);
  return res;
}


int ZGraph::get_latest_not_after_index
(const ZPartialOrder& po, const ZEvent *read, const CPid& cpid) const
{
  assert(is_read(read));
  assert(has_event(read));
  assert(has_thread(cpid));
  assert(cpid != read->cpid());

  assert(_cache.writes.count(read->ml()));
  if (!_cache.writes.at(read->ml()).count(cpid)) // using 'at' because this method is const
    return -1;
  assert(_cache.writes.at(read->ml()).count(cpid));

  int su = po.succ(read, cpid).second;
  if (su >= (int) (*this)(cpid).size()) {
    assert(su == INT_MAX);
    su = (int) (*this)(cpid).size();
  }
  // find tail write in cpid, starting from (and including) su-1, for read-ml
  // return its index in cache.writes
  return get_tailw_index(read->ml(), cpid, su - 1);
}


const ZEvent * ZGraph::get_latest_not_after
(const ZPartialOrder& po, const ZEvent *read, const CPid& cpid) const
{
  assert(is_read(read));
  assert(has_event(read));
  assert(has_thread(cpid));
  assert(cpid != read->cpid());

  int idx = get_latest_not_after_index(po, read, cpid);
  if (idx == -1)
    return nullptr;

  assert(idx < (int) _cache.writes.at(read->ml()).at(cpid).size());
  auto res = _cache.writes.at(read->ml()).at(cpid)[idx]; // using 'at' because this method is const
  assert(is_write(res) && same_ml(res, read) &&
         res->cpid() != read->cpid() && res->cpid() == cpid &&
         !po.has_edge(read, res));
  return res;
}


const ZEvent * ZGraph::get_local_write(const ZEvent *read) const
{
  assert(is_read(read));
  assert(has_event(read));
  assert(_cache.local_write.count(read));
  const ZEvent *local = _cache.local_write.at(read); // using 'at' because this method is const
  assert(!local || (is_write(local) && same_ml(read, local)));
  return local;
}


std::list<const ZEvent *> ZGraph::events_to_mutate(const ZAnnotation& annotation) const
{
  auto res = std::list<const ZEvent *>();
  for (unsigned i = 0; i < _lines.size(); ++i) {
    assert(!_lines[i].empty());
    const ZEvent * last_ev = _lines[i][ _lines[i].size() - 1 ];
    assert(last_ev);
    // When the last event of a thread is lock, it suffices to ask whether it's
    // the 'last-annotated' lock, because if it was annotated before, it would
    // not be the last event of the thread - at least its unlock would follow
    if ((is_read(last_ev) && !annotation.defines(last_ev)) ||
        (is_lock(last_ev) && !annotation.is_last_lock(last_ev)))
      res.push_back(last_ev);
  }
  // Sort the events based on a specified order
  res.sort(ZEventPtrComp());
  return res;
}


std::set<const ZEvent *> ZGraph::mutation_candidates_collect
(const ZPartialOrder& po, const ZEvent *read,
 const std::set<ZEventID>& check_if_any_is_visible) const
{
  assert(is_read(read));
  assert(has_event(read));

  std::set<const ZEvent *> obs_events;

  std::unordered_set<const ZEvent *> notCovered;
  std::unordered_set<const ZEvent *> mayBeCovered;
  // From the value (evid) and below, everything is covered from read by some other write
  std::unordered_map<CPid, int> covered;
  for (const CPid& cpid : threads())
    covered.emplace(cpid, -1);

  // Handle other threads
  for (const CPid& cpid : threads()) {
    if (read->cpid() != cpid) {
      int su = po.succ(read, cpid).second;
      if (su >= (int) (*this)(cpid).size()) {
        assert(su == INT_MAX);
        su = (int) (*this)(cpid).size();
      }
      // Cache-using implementation below; Naive -- see comment after the method
      // Tail write index to cache.writes
      int w_idx = get_tailw_index(read->ml(), cpid, su - 1);
      assert(w_idx == get_latest_not_after_index(po, read, cpid));
      while (w_idx > -1) {
        const ZEvent *rem = _cache.writes.at(read->ml()).at(cpid).at(w_idx); // using 'at' because this method is const
        assert(rem && is_write(rem) && same_ml(read, rem) && rem->cpid() == cpid);
        assert(!po.has_edge(read, rem));
        if (po.has_edge(rem, read)) {
          // Others in this thread are covered from read by rem
          assert(rem->event_id() <= po.pred(read, cpid).second);
          // rem may be covered by some conflicting write rem -> x -> read
          mayBeCovered.emplace(rem);
          // Update cover
          for (const CPid& covcpid : threads()) {
            if (covcpid != rem->cpid()) {
              int newcov = po.pred(rem, covcpid).second;
              assert(covered.count(covcpid));
              if (newcov > covered[covcpid])
                covered[covcpid] = newcov;
            }
          }
          // Break -- Others in this thread are covered from read by rem
          break;
        } else {
          notCovered.emplace(rem);
          --w_idx;
          // VISIBLE check
          if (!check_if_any_is_visible.empty() &&
              check_if_any_is_visible.count(rem->id())) {
            std::set<const ZEvent *> visible_hit;
            visible_hit.emplace(rem);
            return visible_hit;
          }
        }
      }
    }
  }

  // Handle thread of read
  const ZEvent * local = _cache.local_write.at(read);
  if (local) {
    // There is a local write for read
    assert(is_write(local) && same_ml(local, read) &&
           local->cpid() == read->cpid() &&
           local->event_id() < read->event_id());
    // Update cover caused by local
    for (const CPid& covcpid : threads()) {
      if (covcpid != local->cpid()) {
        int newcov = po.pred(local, covcpid).second;
        assert(covered.count(covcpid));
        if (newcov > covered[covcpid])
          covered[covcpid] = newcov;
      }
    }
    // Check if not covered by some remote
    bool localCovered = false;
    for (const auto& rem : mayBeCovered) {
      assert(is_write(rem) && same_ml(rem, read) && po.has_edge(rem, read));
      if (po.has_edge(local, rem)) {
        localCovered = true;
        break;
      }
    }
    if (!localCovered) {
      // Local write not covered, add obs
      obs_events.emplace(local);
      // VISIBLE check
      if (!check_if_any_is_visible.empty() &&
          check_if_any_is_visible.count(local->id())) {
        std::set<const ZEvent *> visible_hit;
        visible_hit.emplace(local);
        return visible_hit;
      }
    }
  } else {
    // No local write for read
    if (mayBeCovered.empty()) {
      // Consider initial event
      assert(initial()->value() == 0);
      obs_events.emplace(initial());
      // VISIBLE check
      if (!check_if_any_is_visible.empty() &&
          check_if_any_is_visible.count(initial()->id())) {
        std::set<const ZEvent *> visible_hit;
        visible_hit.emplace(initial());
        return visible_hit;
      }
    }
  }

  // Take candidates unordered with read
  for (const auto& rem : notCovered) {
    assert(is_write(rem) && same_ml(rem, read) &&
           !po.are_ordered(rem, read));
    obs_events.emplace(rem);
    // VISIBLE check has already been done above for these
    assert(check_if_any_is_visible.empty() ||
           !check_if_any_is_visible.count(rem->id()));
  }

  // Take candidates that happen before read
  for (const auto& rem : mayBeCovered) {
    assert(is_write(rem) && same_ml(rem, read) &&
           po.has_edge(rem, read));
    assert(covered.count(rem->cpid()));
    if (rem->event_id() > covered[rem->cpid()]) {
      // Not covered, add
      obs_events.emplace(rem);
      // VISIBLE check
      if (!check_if_any_is_visible.empty() &&
          check_if_any_is_visible.count(rem->id())) {
        std::set<const ZEvent *> visible_hit;
        visible_hit.emplace(rem);
        return visible_hit;
      }
    }
  }

  // VISIBLE failed if it reached here
  if (!check_if_any_is_visible.empty())
    return std::set<const ZEvent *>();

  return obs_events;
}

      /*
      int curr = su - 1;
      while (curr >= 0) {
        const ZEvent *rem = (*this)(tid, auxid)[curr];
        assert(is_writeM(rem) && !po.has_edge(read, rem));
        if (po.has_edge(rem, read)) {
          if (same_ml(read, rem)) {
            mayBeCovered.emplace(rem);
            // Update cover
            for (unsigned covthr = 0; covthr < number_of_threads(); ++covthr) {
              if (covthr != rem->thread_id()) {
                for (int covaux : auxes(covthr)) {
                  int newcov = po.pred(rem, covthr, covaux).second;
                  auto covta = std::pair<unsigned, int>(covthr, covaux);
                  assert(covered.count(covta));
                  if (newcov > covered[covta])
                    covered[covta] = newcov;
                }
              }
            }
            // Break
            break;
          }
          curr--;
        } else {
          if (same_ml(read, rem))
            notCovered.emplace(rem);
          curr--;
        }
      }
      */

void ZGraph::mutation_candidates_filter_by_negative
(const ZEvent *read, std::set<const ZEvent *>& candidates,
 const ZAnnotationNeg& negative) const
{
  auto it = candidates.begin();
  while (it != candidates.end()) {
    const ZEvent * cand = *it;
    if (is_initial(cand)) {
      if (negative.forbids_initial(read))
        it = candidates.erase(it);
      else
        ++it;
      continue;
    }
    assert(is_write(cand));
    if (negative.forbids(read, cand))
      it = candidates.erase(it);
    else
      ++it;
  }
}


std::set<ZAnn> ZGraph::mutation_candidates_grouped
(const ZEvent *read, const ZAnnotationNeg& negative) const
{
  // Collect the candidates
  std::set<const ZEvent *> obs_events = mutation_candidates_collect
  (_po, read, std::set<ZEventID>());

  // Filter by negative annotation
  mutation_candidates_filter_by_negative(read, obs_events, negative);

  // Group by value
  std::map<int, std::set<ZEventID>> anns;
  for (const ZEvent * ev : obs_events) {
    if (!anns.count(ev->value()))
      anns.emplace(ev->value(), std::set<ZEventID>());
    anns[ev->value()].emplace(ev->id());
  }

  // Create ZAnn-s
  std::set<ZAnn> res;
  for (const auto& v_e : anns)
    res.emplace(ZAnn(v_e.first, v_e.second));
  return res;
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZGraph& gr)
{
  out << gr.to_string();
  return out;
}
