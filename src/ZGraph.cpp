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


// Empty/Initial
ZGraph::ZGraph()
  : _init(ZEvent(true)),
    _lines(),
    _cpid_to_line(),
    _line_to_cpid(),
    _cache()
{
  assert(empty());
}


// Copying/Extending
ZGraph::ZGraph(const ZGraph& oth)
  : _init(ZEvent(true)),
    _lines(oth._lines),
    _cpid_to_line(oth._cpid_to_line),
    _line_to_cpid(oth._line_to_cpid),
    _cache()
{
  assert(!empty());
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


void ZGraph::add_line(const CPid& cpid)
{
  assert(!has_thread(cpid));
  _cpid_to_line.emplace(cpid, _lines.size());
  _line_to_cpid.push_back(cpid);
  assert(_cpid_to_line.size() == _line_to_cpid.size());
  _lines.push_back(LineT());
  _lines.back().reserve(16);
}


void ZGraph::add_line(const ZEvent *ev)
{
  assert(ev && !is_initial(ev));
  assert(!has_event(ev));
  add_line(ev->cpid());
}


void ZGraph::inherit_lines(const ZPartialOrder& po)
{
  assert(empty());
  for (const CPid& cpid : po.threads_spanned()) {
    add_line(cpid);
  }
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


bool ZGraph::has_thread(const CPid& cpid) const
{
  return _cpid_to_line.count(cpid);
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
  int thr_size = po.thread_size(cpid);
  if (su >= thr_size) {
    assert(su == INT_MAX);
    su = thr_size;
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
  assert(po.spans_event(res));
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


std::list<const ZEvent *> ZGraph::events_to_mutate
(const ZPartialOrder& po, const ZAnnotation& annotation) const
{
  auto res = std::list<const ZEvent *>();
  for (unsigned i = 0; i < po.threads_spanned().size(); ++i) {
    assert(i < _lines.size() && i < _line_to_cpid.size());
    const CPid& cpid = line_id_to_cpid(i);
    assert(po.spans_thread(cpid));
    assert(!_lines[i].empty());
    int thr_size = po.thread_size(cpid);
    assert(thr_size > 0);
    const ZEvent * last_ev = _lines[i][ thr_size - 1 ];
    assert(last_ev);
    assert(po.spans_event(last_ev));
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
 const std::set<ZEventID>& check_if_any_is_visible,
 const ZPartialOrder * stricter_po) const
{
  assert(is_read(read));
  assert(has_event(read));
  assert(po.spans_event(read));
  assert(!stricter_po || !stricter_po->spans_event(read));

  std::set<const ZEvent *> obs_events;

  std::unordered_set<const ZEvent *> notCovered;
  std::unordered_set<const ZEvent *> mayBeCovered;
  // From the value (evid) and below, everything is covered from read by some other write
  std::unordered_map<CPid, int> covered;
  for (const CPid& cpid : po.threads_spanned())
    covered.emplace(cpid, -1);

  // Handle other threads
  for (const CPid& cpid : po.threads_spanned()) {
    if (read->cpid() == cpid)
      continue;
    int su = po.succ(read, cpid).second;
    int thr_size = po.thread_size(cpid);
    if (su >= thr_size) {
      assert(su == INT_MAX);
      su = thr_size;
    }
    // Cache-using implementation below; Naive -- see comment after the method
    // Tail write index to cache.writes
    int w_idx = get_tailw_index(read->ml(), cpid, su - 1);
    assert(w_idx == get_latest_not_after_index(po, read, cpid));
    while (w_idx > -1) {
      const ZEvent *rem = _cache.writes.at(read->ml()).at(cpid).at(w_idx); // using 'at' because this method is const
      assert(rem && is_write(rem) && same_ml(read, rem) && rem->cpid() == cpid);
      assert(po.spans_event(rem));
      assert(!po.has_edge(read, rem));
      if (po.has_edge(rem, read)) {
        // Others in this thread are covered from read by rem
        assert(rem->event_id() <= po.pred(read, cpid).second);
        // rem may be covered by some conflicting write rem -> x -> read
        mayBeCovered.emplace(rem);
        // Update cover
        for (const CPid& covcpid : po.threads_spanned()) {
          if (covcpid == rem->cpid())
            continue;
          int newcov = po.pred(rem, covcpid).second;
          assert(covered.count(covcpid));
          if (newcov > covered[covcpid])
            covered[covcpid] = newcov;
        }
        // Break -- Others in this thread are covered from read by rem
        break;
      }
      assert(!po.are_ordered(read, rem));
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

  // Handle thread of read
  const ZEvent * local = _cache.local_write.at(read);
  if (local) {
    // There is a local write for read
    assert(is_write(local) && same_ml(local, read) &&
           local->cpid() == read->cpid() &&
           local->event_id() < read->event_id());
    // Update cover caused by local
    for (const CPid& covcpid : po.threads_spanned()) {
      if (covcpid == local->cpid())
        continue;
      int newcov = po.pred(local, covcpid).second;
      assert(covered.count(covcpid));
      if (newcov > covered[covcpid])
        covered[covcpid] = newcov;
    }
    // Check if not covered by some remote
    assert(covered.count(local->cpid()));
    bool localCovered = (covered.at(local->cpid()) >= local->event_id());
  #ifndef NDEBUG
    bool loc = false;
    for (const auto& rem : mayBeCovered) {
      assert(is_write(rem) && same_ml(rem, read) && po.has_edge(rem, read));
      if (po.has_edge(local, rem)) {
        loc = true;
        break;
      }
    }
    assert((!localCovered || loc) && (localCovered || !loc));
  #endif
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


void ZGraph::mutation_candidates_filter_by_stricter
(std::set<const ZEvent *>& candidates,
 const ZPartialOrder& stricter_po) const
{
  auto it = candidates.begin();
  while (it != candidates.end()) {
    const ZEvent * cand = *it;
    if (is_initial(cand)) {
      it = candidates.erase(it);
      continue;
    }
    assert(is_write(cand));
    if (stricter_po.spans_event(cand))
      it = candidates.erase(it);
    else
      ++it;
  }
}


std::set<ZAnn> ZGraph::mutation_candidates_grouped
(const ZPartialOrder& po, const ZEvent *read,
 const ZAnnotationNeg& negative,
 const ZPartialOrder * stricter_po) const
{
  // Do not use stricter_po to collect, in case it spans read
  bool stricter_spans_read = (stricter_po && stricter_po->spans_event(read));
  const ZPartialOrder * stricter_to_collect = (
  stricter_spans_read ? nullptr : stricter_po);

  // Collect the candidates
  std::set<const ZEvent *> obs_events = mutation_candidates_collect
  (po, read, std::set<ZEventID>(), stricter_to_collect);

  // Filter by negative annotation
  mutation_candidates_filter_by_negative(read, obs_events, negative);

  // If stricter_spans_read, delete all events spanned by stricter_po
  if (stricter_spans_read) {
    mutation_candidates_filter_by_stricter(obs_events, *stricter_po);
  }

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
