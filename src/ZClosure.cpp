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

#include "ZClosure.h"


/* *************************** */
/* RULE 2                      */
/* *************************** */

// if badwrite -> read, then badwrite -> write
// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::rule_two
(const ZEvent *read, const ZEvent *write)
{
  assert(read && is_read(read));
  assert(write && (is_initial(write) ||
                   (is_write(write) && same_ml(read, write))));
  assert(!po.has_edge(read, write));

  bool change = false;
  // Thread of read was handled in rule 1
  // Handle thread of write
  if (!is_initial(write) && write->cpid() != read->cpid()) {
    int last_pred = po.pred(read, write->cpid()).second;
    if (last_pred >= 0) {
      const ZEvent *last_conf_pred = gr.get_tailw(read->ml(),
                                                  write->cpid(), last_pred);
      if (last_conf_pred) {
        assert(last_conf_pred && same_ml(read, last_conf_pred) &&
               po.has_edge(last_conf_pred, read));
        assert(write->cpid() == last_conf_pred->cpid());
        if (write->event_id() < last_conf_pred->event_id())
          return {true, false}; // Impossible - covered by same-thread event
      }
    }
  }

  // Handle all threads except thread of read and thread of write
  for (const CPid& cpid : po.threads_spanned()) {
    if (cpid == read->cpid() || (!is_initial(write) && cpid == write->cpid()))
      continue;
    int last_pred = po.pred(read, cpid).second;
    if (last_pred == -1)
      continue;
    const ZEvent *last_conf_pred = gr.get_tailw(read->ml(), cpid, last_pred);
    if (!last_conf_pred)
      continue;
    assert(last_conf_pred && same_ml(read, last_conf_pred) &&
           po.has_edge(last_conf_pred, read));
    if (is_initial(write) || po.has_edge(write, last_conf_pred))
      return {true, false}; // Impossible - reverse edge already present
    if (!po.has_edge(last_conf_pred, write)) {
      po.add_edge(last_conf_pred, write);
      change = true;
      added_edges++;
    }
  }

  return {false, change}; // done, change-variable
}

/* *************************** */
/* RULE 3                      */
/* *************************** */

// if write -> badwrite, then read -> badwrite
// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::rule_three
(const ZEvent *read, const ZEvent *write) {
  assert(read && is_read(read));
  assert(write && (is_initial(write) ||
                   (is_write(write) && same_ml(read, write))));
  assert(!po.has_edge(read, write));

  bool change = false;
  // Handle thread of read
  const ZEvent * local = gr.get_local_write(read);
  assert(!local || po.spans_event(local));
  if (local && (*local) != (*write) && po.has_edge(write, local))
    return {true, false}; // Impossible - covered by read-local badwrite

  // Handle thread of write
  if (!is_initial(write) && write->cpid() != read->cpid()) {
    int idx = gr.get_latest_not_after_index(po, read, write->cpid());
    // {0,1,..,idx} in cache at same ml not happening after read
    if (idx != -1) {
      const auto& wrcache = gr.cache().writes.at(read->ml()).at(write->cpid());
      assert(idx < (int) wrcache.size());

      // Binary search in cache events to get first badwrite after write
      int l = 0;
      int r = idx;
      while(l<r) {
        int mid=(l+r)>>1;
        auto res = wrcache.at(mid);
        assert(is_write(res) && same_ml(res, read) &&
               res->cpid() == write->cpid());
        if (res->event_id() > write->event_id()) r=mid;
        else l=mid+1;
      } // after the loop, l = r
      assert(l == r);
      auto res = wrcache.at(l);
      assert(is_write(res) && same_ml(res, read) &&
            res->cpid() == write->cpid());
      if (res->event_id() > write->event_id()) {
        if (po.has_edge(res, read))
          return {true, false}; // Impossible - reverse edge already present
        if (!po.has_edge(read, res)) {
          po.add_edge(read, res);
          change = true;
          added_edges++;
        }
      }
    }
  }

  // Handle all threads except thread of read and thread of write
  for (const CPid& cpid : po.threads_spanned()) {
    if (cpid == read->cpid() || (!is_initial(write) && cpid == write->cpid()))
      continue;
    int idx = gr.get_latest_not_after_index(po, read, cpid);
    // {0,1,..,idx} in cache at same ml not happening after read
    if (idx == -1)
      continue;
    const auto& wrcache = gr.cache().writes.at(read->ml()).at(cpid);
    assert(idx < (int) wrcache.size());

    // Binary search in cache events to get first badwrite after write
    int l = 0;
    int r = idx;
    while(l<r) {
      int mid=(l+r)>>1;
      auto res = wrcache.at(mid);
      assert(is_write(res) && same_ml(res, read) &&
             res->cpid() == cpid);
      if (po.has_edge(write, res)) r=mid;
      else l=mid+1;
    } // after the loop, l = r
    assert(l == r);
    auto res = wrcache.at(l);
    assert(is_write(res) && same_ml(res, read) &&
           res->cpid() == cpid);
    if (po.has_edge(write, res)) {
      if (po.has_edge(res, read))
        return {true, false}; // Impossible - reverse edge already present
      if (!po.has_edge(read, res)) {
        po.add_edge(read, res);
        change = true;
        added_edges++;
      }
    }
  }

  return {false, change}; // done, change-variable
}

/* *************************** */
/* RULES                       */
/* *************************** */

std::pair<bool, bool> ZClosure::rules
(const ZEvent *read, const ZAnn& ann)
{
  assert(read && is_read(read));
  assert(ann.goodwrites.size() == 1);
  const ZEvent * write = gr.event(*ann.goodwrites.begin());
  assert(po.spans_event(write));
  bool change = false;
  // Rule2
  auto res = rule_two(read, write);
  if (res.first) return {true, false};
  if (res.second) change = true;
  // Rule3
  res = rule_three(read, write);
  if (res.first) return {true, false};
  if (res.second) change = true;
  return {false, change};
}

/* *************************** */
/* CLOSE                       */
/* *************************** */

bool ZClosure::close_finish
(const std::map<const ZEvent *, ZAnn>& reads_with_multiple_good_writes) const
{
  for (const auto& read_ann : reads_with_multiple_good_writes) {
    assert(read_ann.second.goodwrites.size() >= 2);
    assert(read_ann.first && is_read(read_ann.first));
    std::set<const ZEvent *> any_good_write_visible = (
      gr.mutation_candidates_collect
      (po, read_ann.first, read_ann.second.goodwrites, nullptr));
    if (any_good_write_visible.empty())
      return false;
    else {
      assert(any_good_write_visible.size() == 1);
      assert(*(any_good_write_visible.begin()));
      assert(read_ann.second.goodwrites.count(
             (*(any_good_write_visible.begin()))->id()));
    }
  }
  return true;
}


bool ZClosure::close(const ZEvent *newread) {
  // Reads with multiple good writes will be checked at the end
  std::map<const ZEvent *, ZAnn> reads_with_multiple_good_writes;

  // Rules for new read
  if (newread) {
    assert(is_read(newread) && po.spans_event(newread));
    const ZAnn& newread_ann = an.ann(newread);
    // If 1 good write perform rules, otherwise check at the end
    if (newread_ann.goodwrites.size() != 1) {
      assert(newread_ann.goodwrites.size() >= 2);
      reads_with_multiple_good_writes.emplace(newread, newread_ann);
    } else {
      auto res = rules(newread, newread_ann);
      if (res.first) { return false; }
    }
  }

  bool change = true;
  int last_change = an.size();
  while (change) {
    iterations++;
    change = false;
    int cur = -1;
    for (const auto& read_ann : an) {
      cur++;
      if (cur == last_change) {
        // One entire iteration without any changes
        // hence the partial order is closed
        assert(!change);
        return close_finish(reads_with_multiple_good_writes);
      }
      const ZEvent * read = gr.event(read_ann.first);
      assert(po.spans_event(read));
      if (newread && (*newread) == (*read))
        continue;
      if (po.is_closure_safe(read))
        continue;
      // If 1 good write perform rules, otherwise check at the end
      if (read_ann.second.goodwrites.size() != 1) {
        assert(read_ann.second.goodwrites.size() >= 2);
        if (!reads_with_multiple_good_writes.count(read))
          reads_with_multiple_good_writes.emplace(read, read_ann.second);
      } else {
        auto res = rules(read, read_ann.second);
        if (res.first) { return false; }
        if (res.second) {
          change = true;
          last_change = cur;
        }
      }
    }
    cur++;
    assert(cur == (int) an.size());
    if (cur == last_change) {
      // One entire iteration without any changes
      // hence the partial order is closed
      assert(!change);
      return close_finish(reads_with_multiple_good_writes);
    }
    if (newread) {
      assert(is_read(newread));
      const ZAnn& newread_ann = an.ann(newread);
      // If 1 good write perform rules,
      // otherwise already added beforehand to be checked at the end
      if (newread_ann.goodwrites.size() == 1) {
        auto res = rules(newread, an.ann(newread));
        if (res.first) { return false; }
        if (res.second) {
          change = true;
          last_change = cur;
        }
      }
    }
  }
  assert(!change);
  return close_finish(reads_with_multiple_good_writes);
}

/* *************************** */
/* RULE 1                      */
/* *************************** */

void ZClosure::rule_one_lock(const ZEvent *lock, const ZEvent *unlock)
{
  assert(is_lock(lock) && is_unlock(unlock) && same_ml(lock, unlock));
  assert(!po.has_edge(lock, unlock));
  if (!po.has_edge(unlock, lock))
    po.add_edge(unlock, lock);
}

void ZClosure::rule_one_multi_good(const ZEvent *read, const ZAnn& ann)
{
  assert(is_read(read) && po.spans_event(read));
  assert(ann.goodwrites.size() >= 2);
  // If there is a po-smallest write, it should happen before read
  std::set<const ZEvent *> good;
  const ZEvent * candidate = nullptr;
  for (const auto& goodwrite_id : ann.goodwrites) {
    const ZEvent * write = gr.event(goodwrite_id);
    assert(po.spans_event(write));
    if (is_initial(write))
      return; // Initial is po-smallest and already happens before read
    assert(is_write(write));
    good.emplace(write);
    if (!candidate || po.has_edge(write, candidate))
      candidate = write;
  }
  // candidate holds a po-minimal write, check if it is po-smallest
  assert(candidate);
  for (const auto& goodwrite : good) {
    if (goodwrite == candidate)
      continue;
    if (!po.has_edge(candidate, goodwrite))
      return; // There is no po-smallest write
  }
  // candidate holds a po-smallest write
  #ifndef NDEBUG
  for (const auto& goodwrite_id : ann.goodwrites) {
    const ZEvent * write = gr.event(goodwrite_id);
    assert(candidate == write || po.has_edge(candidate, write));
  }
  #endif
  if (!po.has_edge(candidate, read))
    po.add_edge(candidate, read);
}

void ZClosure::rule_one(const ZEvent *read, const ZAnn& ann)
{
  assert(is_read(read) && po.spans_event(read));
  if (ann.goodwrites.size() != 1) {
    assert(ann.goodwrites.size() >= 2);
    rule_one_multi_good(read, ann);
  }
  const ZEvent * write = gr.event(*ann.goodwrites.begin());
  assert(po.spans_event(write));
  if (is_initial(write)) {
    assert(!gr.get_local_write(read));
    return;
  }
  assert(is_write(write) && same_ml(read, write));
  if (read->cpid() != write->cpid()) {
    assert(!po.has_edge(read, write));
    if (!po.has_edge(write, read))
      po.add_edge(write, read);
    // Cover local write from read by local -> write -> read
    const ZEvent * local = gr.get_local_write(read);
    if (local) {
      assert(is_write(local) && same_ml(local, read) &&
             local->cpid() == read->cpid() &&
             local->event_id() < read->event_id());
      assert(!po.has_edge(write, local));
      if (!po.has_edge(local, write))
        po.add_edge(local, write);
    }
  } else
    assert(write == gr.get_local_write(read));
}
