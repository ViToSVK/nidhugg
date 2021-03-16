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
#include <iostream>


ZClosure::ZClosure
(const ZAnnotation& annotation, ZGraph& graph, ZPartialOrder& partial)
  : an(annotation),
    gr(graph),
    po(partial) { assert(&(po.graph) == &graph); }


// return writeBuffer, writeMemory for obs
std::pair<const ZEvent *, const ZEvent *> ZClosure::getObs(const ZEventID& id) {
  if(id.event_id() == -1)
    return {nullptr,nullptr};
  auto buffer_part = gr.getEvent(id);
  assert(isWriteB(buffer_part));
  auto write_part = buffer_part->write_other_ptr;
  assert(isWriteM(write_part));
  return {buffer_part,write_part};
}


/* *************************** */
/* RULE 1                      */
/* *************************** */

// true iff impossible
bool ZClosure::ruleOne(const ZEvent *ev, const ZEventID& obs) {
  assert(gr.hasEvent(ev));
  if (obs.event_id() == -1) {
    // Handle initial-event observation separately
    // Nothing to do in rule1
    return false;
  }

  const ZEvent *obsEv = gr.getEvent(obs);
  assert(obsEv && sameMl(ev, obsEv));
  assert((isRead(ev) && isWriteB(obsEv)) ||
         (isLock(ev) && isUnlock(obsEv)));
  assert(gr.hasEvent(obsEv));

  if (isLock(ev)) {
    // Handle lock annotation
    if (po.hasEdge(ev, obsEv))
      return true;
    if (!po.hasEdge(obsEv, ev))
      po.addEdge(obsEv, ev);
    return false;
  }

  // Handle read annotation
  assert(isRead(ev) && obsEv->write_other_ptr);
  if (ev->thread_id() != obsEv->thread_id()) {
    auto obsMem = obsEv->write_other_ptr;
    assert(isWriteM(obsMem));
    assert(gr.hasEvent(obsMem));
    if (po.hasEdge(ev, obsMem))
      return true;
    if (!po.hasEdge(obsMem, ev))
      po.addEdge(obsMem, ev);
    // Edge to 'obsMem' from memory-write of largest (by eventid)
    // local buffer write which happens before read 'ev'
    auto lastBuf = gr.getLocalBufferW(ev);
    if(lastBuf) {
      assert(isWriteB(lastBuf));
      assert(gr.hasEvent(lastBuf));
      auto mem_counterpart = lastBuf->write_other_ptr; // Getting the Memory Write
      assert(isWriteM(mem_counterpart));
      assert(gr.hasEvent(mem_counterpart));
      if (po.hasEdge(obsMem, mem_counterpart))
        return true;
      if(!po.hasEdge(mem_counterpart, obsMem)) {
        po.addEdge(mem_counterpart, obsMem); // Memory write before curr obs
        added_edges++;
      }
    }
  }
  return false;
}


/* *************************** */
/* RULE 2                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleTwo(const ZEvent *read, const ZEventID& obs) {
  // Optimization idea: if obs different thread then only once this function could be called
  // as no memory write will be added in between other rules
  // To do so we can call it preclose and then in rules don't call again

  auto write = getObs(obs);
  auto write_memory = write.second;
  assert(!write_memory || isWriteM(write_memory));
  assert(!write_memory || gr.hasEvent(write_memory));
  // Idea: Iterate on all (but self) threads, in each thread i:
  // get 'memory_pred' - conflicting memory-write predecessor of read in thread i
  // add edge if not already present, that 'memory_pred' -> 'write_memory'
  //unsigned totThreads =  gr.number_of_threads();
  unsigned readThread = read->thread_id();
  unsigned writeThread = write_memory ? write_memory->thread_id() : INT_MAX;
  auto location = read->ml; // SymAddrSize of read
  bool change = false;
  for(unsigned i : gr.get_threads()) { // looping over all threads
    if(i != readThread and i != writeThread) { // all but read thread and observation-write thread
      int readaux = gr.auxForMl(location, i);
      if(readaux == -1) continue;
      int lastEvid = po.pred(read,i,readaux).second ; // Hardcode Aux gets Event id of pred
      if(lastEvid != -1) { // {0,1,2,..,lastEvid} happen before read
        auto memory_pred = gr.getTailW(location,i,lastEvid); // last write of thread i on same location as read
        if(memory_pred) { // not nullptr
          assert(isWriteM(memory_pred));
          assert(gr.hasEvent(memory_pred));
          if(!write_memory || po.hasEdge(write_memory,memory_pred))
            return {true,false}; // Impossible - reverse edge already present
          if(!po.hasEdge(memory_pred,write_memory)) {
            po.addEdge(memory_pred,write_memory);
            change = true;
            added_edges++;
          }
        }
      }
    }
  }
  if(write_memory and readThread != writeThread) { // checking for impossibility due to write thread
    int writeaux = gr.auxForMl(location, writeThread);
    if(writeaux != -1) {
      int lastEvid = po.pred(read,writeThread,writeaux).second ; // Hardcode Aux gets Event id of pred
      if(lastEvid != -1) { // {0,1,2,..,lastEvid} happen before read
        auto memory_pred = gr.getTailW(location,writeThread,lastEvid);  // last write of thread i on same location as read
        if(memory_pred) { // not nullptr
          assert(isWriteM(memory_pred));
          assert(gr.hasEvent(memory_pred));
          if(memory_pred->event_id() > write_memory->event_id())
            return {true,false}; // Impossible - reverse edge already present
        }
      }
    }
  }
  //po.dump();
  return {false, change}; // done, change-variable
}

  // What type of edges are added ?
    // Rule 1 : obs memory on read (diff thread) , self memory to read
    // Rule 2 : memory to obs memory (non read thread) if memory < read
    // Rule 3 : read to memory (non read thread) if obs memory < memory
  // Rule 3 has to be involved again for any type
  // Seems like Rule 2 happens only once


/* *************************** */
/* RULE 3                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleThree(const ZEvent *read, const ZEventID& obs) {
  // Without optimization
  auto write = getObs(obs);
  auto write_memory = write.second;
  assert(!write_memory || isWriteM(write_memory));
  // Idea: Iterate on all (but self and writeM) threads, for each thread i:
  // get 'bad_writeM' - conflicting memory-write successor of write_memory
  // add edge if not already present, that 'read' -> 'bad_writeM'
  //unsigned totThreads =  gr.number_of_threads();
  unsigned readThread = read->thread_id();
  unsigned writeThread = write_memory ? write_memory->thread_id():INT_MAX;
  bool change = false;
  for(unsigned i : gr.get_threads()) {  // looping over all threads
    if(i != readThread and i != writeThread) {  // all but read thread
      int lastBefore = gr.getLatestNotAfterIndex(read,i,po); // {0,1,..,lastBefore} in cache at same ml with r </ wM
      if(lastBefore == -1) continue;
      auto &cache = gr.getCache();
      assert(lastBefore < (int) cache.wm.at(read->ml).at(i).size());
      if(write_memory) {
        int l=0,r=lastBefore;
        // Binary search in cache events to get first wM after write
        while(l<r) { // gives r if fail
          int mid=(l+r)>>1;
          auto res = cache.wm.at(read->ml).at(i)[mid];
          assert(isWriteM(res) && sameMl(res, read) &&
                 res->thread_id() == i && res->aux_id() != -1);
          if(po.hasEdge(write_memory,res)) r=mid;
          else l=mid+1;
        } // after the loop, l = r = x
        auto res = cache.wm.at(read->ml).at(i)[l];
        assert(isWriteM(res) && sameMl(res, read) &&
               res->thread_id() == i && res->aux_id() != -1);
        if(po.hasEdge(write_memory,res)) {
          if(po.hasEdge(res,read))
            return {true,false}; // Impossible - reverse edge already present
          if(!po.hasEdge(read,res)) {
            po.addEdge(read,res);
            change = true;
            added_edges++;
          }
        }
      } else {  // initial-event observation
        auto res = cache.wm.at(read->ml).at(i)[0]; // adding edge from first one
        assert(isWriteM(res) && sameMl(res, read) &&
               res->thread_id() == i && res->aux_id() != -1);
        if(po.hasEdge(res,read))
          return {true,false}; // Impossible - reverse edge already present
        if(!po.hasEdge(read,res)) {
          po.addEdge(read,res);
          change = true;
          added_edges++;
        }
      }
    }
  }
  // Proceed as above, on the thread of 'write_memory'
  if(write_memory and writeThread != readThread) {
    int lastBefore = gr.getLatestNotAfterIndex(read,writeThread,po); // {0,1,..,lastBefore} in cache at same ml with r </ wM
    if(lastBefore != -1) {
      auto cache = gr.getCache();
      assert(lastBefore < (int) cache.wm.at(read->ml).at(writeThread).size());
      // Binary search in cache events to get first wM after write
      int l=0,r=lastBefore;
      unsigned memId = write_memory->event_id();
      while(l<r) { // gives r if fail
        int mid=(l+r)>>1;
        auto res = cache.wm.at(read->ml).at(writeThread)[mid];
        assert(isWriteM(res) && sameMl(res, read) &&
            res->thread_id() == writeThread && res->aux_id() != -1);
        if(res->event_id()>memId) r=mid;
        else l=mid+1;
      } // after the loop, l = r = x
      auto res = cache.wm.at(read->ml).at(writeThread)[l];
      assert(isWriteM(res) && sameMl(res, read) &&
            res->thread_id() == writeThread && res->aux_id() != -1);
      if(res->event_id()>memId) {
        if(po.hasEdge(res,read))
          return {true,false}; // Impossible - reverse edge already present
        if(!po.hasEdge(read,res)) {
          po.addEdge(read,res);
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

std::pair<bool, bool> ZClosure::rulesTwoThree
(const ZEvent *read, const ZEventID& obs) {
  assert(read && isRead(read));
  assert(gr.hasEvent(read));
  bool change = false;
  // Rule1 is done only the first time
  // Rule2
  auto res = ruleTwo(read, obs);
  if (res.first) return {true, false};
  if (res.second) change = true;
  //Rule3
  res = ruleThree(read, obs);
  if (res.first) return {true, false};
  if (res.second) change = true;
  return {false, change};
}


/* *************************** */
/* CLOSE                       */
/* *************************** */
bool ZClosure::close() {
  // RULE 1 FOR ALL READS
  for (auto it = an.read_begin(); it != an.read_end(); ++it) {
    auto read = gr.getEvent(it->first);
    assert(isRead(read));
    bool res = ruleOne(read, it->second);
    if (res) { return false; }
  }
  // RULE 1 FOR ALL LOCKS
  for (auto it = an.lock_begin(); it != an.lock_end(); ++it) {
    auto lock = gr.getEvent(it->first);
    assert(isLock(lock));
    bool res = ruleOne(lock, it->second);
    if (res) { return false; }
  }

  // IN A LOOP, RULE 2+3 FOR ALL READS UNTIL NO CHANGE
  bool change = true;
  int last_change = an.read_size();
  while (change) {
    iterations++;
    change = false;
    int cur = -1;
    for (auto it = an.read_begin(); it != an.read_end(); ++it) {
      cur++;
      if (cur == last_change) {
        // One entire iteration without any changes
        // hence the partial order is closed
        return true;
      }
      auto read = gr.getEvent(it->first);
      assert(isRead(read));
      auto res = rulesTwoThree(read, it->second);
      if (res.first) { return false; }
      if (res.second) {
        change = true;
        last_change = cur;
      }
    }
    cur++;
    assert(cur == (int) an.read_size());
    if (cur == last_change) {
      // One entire iteration without any changes
      // hence the partial order is closed
      return true;
    }
  }
  assert(!change);
  return true;
}
