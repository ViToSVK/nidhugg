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
  for (const CPid& cpid : gr.threads()) {
    if (cpid != read->cpid() && (is_initial(write) || cpid != write->cpid())) {
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
  }

  //po.dump();
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
  if (local && (*local) != (*write) && po.has_edge(write, local))
    return {true, false}; // Impossible - covered by read-local badwrite

  // Handle thread of write
  if (!is_initial(write) && write->cpid() != read->cpid()) {
    int idx = gr.get_latest_not_after_index(read, write->cpid());
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
          return {true,false}; // Impossible - reverse edge already present
        if (!po.has_edge(read, res)) {
          po.add_edge(read, res);
          change = true;
          added_edges++;
        }
      }
    }
  }

  // Handle all threads except thread of read and thread of write
  for (const CPid& cpid : gr.threads()) {
    if (cpid != read->cpid() && (is_initial(write) || cpid != write->cpid())) {
      int idx = gr.get_latest_not_after_index(read, cpid);
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
  }

  // po.dump();
  return {false, change}; // done, change-variable
}

/* *************************** */
/* RULES                       */
/* *************************** */

std::pair<bool, bool> ZClosure::rules
(const ZEvent *read, const ZAnn& ann)
{
  assert(read && is_read(read));
  if (ann.goodwrites.size() != 1) {
    assert(ann.goodwrites.size() >= 2);
    return {false, false};
  }
  const ZEvent * write = gr.event(*ann.goodwrites.begin());
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

bool ZClosure::close(const ZEvent *newread) {
  // Rules for new read
  if (newread) {
    assert(is_read(newread));
    auto res = rules(newread, an.ann(newread));
    if (res.first) { return false; }
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
        return true;
      }
      const ZEvent * read = gr.event(read_ann.first);
      if (newread && (*newread) == (*read))
        continue;
      if (po.is_closure_safe(read))
        continue;
      auto res = rules(read, read_ann.second);
      if (res.first) { return false; }
      if (res.second) {
        change = true;
        last_change = cur;
      }
    }
    cur++;
    assert(cur == (int) an.size());
    if (cur == last_change) {
      // One entire iteration without any changes
      // hence the partial order is closed
      assert(!change);
      return true;
    }
    if (newread) {
      auto res = rules(newread, an.ann(newread));
      if (res.first) { return false; }
      if (res.second) {
        change = true;
        last_change = cur;
      }
    }
  }
  //po.dump();
  //an.dump();
  assert(!change);
  return true;
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

void ZClosure::rule_one(const ZEvent *read, const ZAnn& ann)
{
  assert(is_read(read));
  if (ann.goodwrites.size() != 1) {
    assert(ann.goodwrites.size() >= 2);
    return;
  }
  const ZEvent * write = gr.event(*ann.goodwrites.begin());
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
