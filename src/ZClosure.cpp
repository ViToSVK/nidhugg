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

#include "ZClosure.h"
#include <iostream>
// Got ZHelper due to header -> Graph -> AnnotationNeg -> Annotation -> Helper

// return writeBuffer, writeMemory for obs
std::pair<const ZEvent *, const ZEvent *> ZClosure::getObs(const ZObs& obs){
	if(obs.thr == INT_MAX) return {nullptr,nullptr};
	auto buffer_part = ba.getEvent(obs.thr, -1, obs.ev);
	assert(isWriteB(buffer_part));
	auto write_part = buffer_part->write_other_ptr;
	assert(isWriteM(write_part));
	return {buffer_part,write_part};
}


/* *************************** */
/* RULE 1                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleOne(const ZEvent *read, const ZObs& obs){
	// Seems like done in preclose
	// Recheck of impossible or change
  return {false, false}; // done, no change
}


/* *************************** */
/* RULE 2                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
//# warning Initial event is not special case right?
std::pair<bool, bool> ZClosure::ruleTwo(const ZEvent *read, const ZObs& obs){
	// Optimization idea if obs different thread then only once this function could be called
	// as no memory write will be added in between other rules
	// To do so we can call it preclose and then in rules don't call again

	// Without optimization
	auto write = getObs(obs);
	auto write_memory = write.second;
	assert(!write_memory || isWriteM(write_memory));
	// Idea: Iterate on all (but self) threads and get memory pred of r
	// add cond that writeM occur after aboce memory write
	unsigned totThreads =  ba.number_of_threads();
	unsigned readThread = read->threadID();
	unsigned writeThread = write_memory ? write_memory->threadID():INT_MAX;
	auto location = read->ml; // SymAddrSize of read
	bool change = false;
	for(unsigned i=0; i < totThreads; ++i){  // looping over all threads
		if(i != readThread and i != writeThread){	// all but read thread
			int readaux = ba.auxForMl(location, i);
			if(readaux == -1) continue;
			int lastEvid = po.pred(read,i,readaux).second ; // Hardcode Aux gets Event id of pred
			if(lastEvid != -1){ // {0,1,2,..,lastEvid} happens before real
				auto memory_pred = gr.getTailW(location,i,lastEvid);	// last write of thread i on same location as read
				if(memory_pred){	// not nullptr
					assert(isWriteM(memory_pred));
					if(!write_memory || po.hasEdge(write_memory,memory_pred)) return {true,false}; // Already reverse edge
        			if(!po.hasEdge(memory_pred,write_memory)){
						po.addEdge(memory_pred,write_memory);
						change = true;
						added_edges++;
					}
				}
			}
		}
	}
	if(write_memory and readThread != writeThread){ // checking for impossibility due to write thread
		int writeaux = ba.auxForMl(location, writeThread);
		if(writeaux != -1){
			int lastEvid = po.pred(read,writeThread,writeaux).second ; // Hardcode Aux gets Event id of pred
			if(lastEvid != -1){ // {0,1,2,..,lastEvid} happens before real
				auto memory_pred = gr.getTailW(location,writeThread,lastEvid);	// last write of thread i on same location as read
				if(memory_pred){	// not nullptr
					assert(isWriteM(memory_pred));
					if(memory_pred->eventID() > write_memory->eventID()) return {true,false}; // Already reverse edge
				}
			}
		}
	}
	//std::cout<<"Rule2 "<<change<<"\n";
	//po.dump();
	return {false, change}; // done, no change
}

	// What type of edges are added ?
		// Rule 1 : obs memory on read (diff thread) , self memory to read
		// Rule 2 : memory to obs memory (non read thread) if memory < read
		// Rule 3 : read to memory (non read thread) if obs memory < memory
	// Rule 3 has to involved again for any type
	// Seems like Rule 2 happens only once


/* *************************** */
/* RULE 3                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> ZClosure::ruleThree(const ZEvent *read, const ZObs& obs){
	// Without optimization
	auto write = getObs(obs);
	auto write_memory = write.second;
	assert(!write_memory || isWriteM(write_memory));
	// Idea: Iterate on all (but self and writeM) threads and get memory suc of writeM
	// add cond that read occur before above memory write
	// Do same stuff for first suc of writeM using different way
	unsigned totThreads =  ba.number_of_threads();
	unsigned readThread = read->threadID();
	unsigned writeThread = write_memory ? write_memory->threadID():INT_MAX;
	bool change = false;
	for(unsigned i=0; i < totThreads; ++i){  // looping over all threads
		if(i != readThread and i != writeThread){	// all but read thread
			int lastBefore = gr.getLatestNotAfterIndex(read,i,po);	// {0,1,..,lastBefore} in cache at same ml with r </ wM
			if(lastBefore == -1) continue;
			auto &cache = gr.getCache();
			assert(lastBefore < (int) cache.wm.at(read->ml).at(i).size());
			if(write_memory){
				int l=0,r=lastBefore;
				// Binary search in cache events to get first wM after write
				while(l<r){ // gives r if fail
					int mid=(l+r)>>1;
					auto res = cache.wm.at(read->ml).at(i)[mid];
					assert(isWriteM(res) && sameMl(res, read) &&
					    res->threadID() == i && res->auxID() != -1);
					if(po.hasEdge(write_memory,res)) r=mid;
					else l=mid+1;
				} // after the loop, l = r = x
				auto res = cache.wm.at(read->ml).at(i)[l];
				assert(isWriteM(res) && sameMl(res, read) &&
					    res->threadID() == i && res->auxID() != -1);
				//res->dump();
				//std::cout<<i<<" "<<std::endl;
				if(po.hasEdge(write_memory,res)){
					if(po.hasEdge(res,read)) return {true,false}; // Already reverse edge
					if(!po.hasEdge(read,res)){
						po.addEdge(read,res);
						change = true;
						added_edges++;
					}
				}
			}else{	// initial obs
				auto res = cache.wm.at(read->ml).at(i)[0]; // adding edge from first one
				assert(isWriteM(res) && sameMl(res, read) &&
					    res->threadID() == i && res->auxID() != -1);
				if(po.hasEdge(res,read)) return {true,false}; // Already reverse edge
				if(!po.hasEdge(read,res)){
					po.addEdge(read,res);
					change = true;
					added_edges++;
				}
			}
		}
	}
	if(write_memory and writeThread != readThread){
		int lastBefore = gr.getLatestNotAfterIndex(read,writeThread,po);// {0,1,..,lastBefore} in cache at same ml with r </ wM
		if(lastBefore != -1){
			auto cache = gr.getCache();
			assert(lastBefore < (int) cache.wm.at(read->ml).at(writeThread).size());
			// Binary search in cache events to get first wM after write
			int l=0,r=lastBefore;
			unsigned memId = write_memory->eventID();
			while(l<r){ // gives r if fail
				int mid=(l+r)>>1;
				auto res = cache.wm.at(read->ml).at(writeThread)[mid];
				assert(isWriteM(res) && sameMl(res, read) &&
				    res->threadID() == writeThread && res->auxID() != -1);
				if(res->eventID()>memId) r=mid;
				else l=mid+1;
			} // after the loop, l = r = x
			auto res = cache.wm.at(read->ml).at(writeThread)[l];
			assert(isWriteM(res) && sameMl(res, read) &&
				    res->threadID() == writeThread && res->auxID() != -1);
			if(res->eventID()>memId){
				if(po.hasEdge(res,read)) return {true,false}; // Already reverse edge
				if(!po.hasEdge(read,res)){
					po.addEdge(read,res);
					change = true;
				}
			}
		}
	}
	//std::cout<<"Rule3 "<<change<<"\n";
	// po.dump();
	return {false, change}; // done, no change
}


/* *************************** */
/* RULES                       */
/* *************************** */

std::pair<bool, bool> ZClosure::rules(const ZEvent *read, const ZObs& obs){
  assert(read && isRead(read));
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
bool ZClosure::close(const ZEvent *newread){
  // Rules for new read
  if (newread) {
    auto res = rules(newread, an.getObs(newread));
    if (res.first) { return false; }
  }

  bool change = true;
  while (change) {
    iterations++;
    change = false;
    for (const auto& read_obs : an) {
      auto read = ba.getEvent(read_obs.first.thr, -1, read_obs.first.ev);
      if (!newread || newread != read) {
        auto res = rules(read, read_obs.second);
        if (res.first) { return false; }
        if (res.second) change = true;
      }
    }
    if (newread) {
      auto res = rules(newread, an.getObs(newread));
      if (res.first) { return false; }
      if (res.second) change = true;
    }
  }
  //po.dump();
  //an.dump();
  assert(!change);
  return true;
}


/* *************************** */
/* PRE-CLOSE                   */
/* *************************** */

void ZClosure::preClose(const ZEvent *ev, const ZEvent *obsEv){
  assert(sameMl(ev, obsEv));
  assert((isRead(ev) && isWriteB(obsEv)) ||
         (isLock(ev) && isUnlock(obsEv)));

  if (isLock(ev)) {
    assert(!po.hasEdge(ev, obsEv));
    if (!po.hasEdge(obsEv, ev)) po.addEdge(obsEv, ev);
    return;
  }

  assert(isRead(ev) && obsEv->write_other_ptr);
  if (ev->threadID() != obsEv->threadID()) {
    auto obsMem = obsEv->write_other_ptr;
    assert(isWriteM(obsMem));
    assert(!po.hasEdge(ev, obsMem));
    if (!po.hasEdge(obsMem, ev)) po.addEdge(obsMem, ev);
	// Adding edge from Memory write of largest (by eventid) Buffer write which happended before current read i.e. ev
    auto lastBuf = gr.getLocalBufferW(ev);
    if(lastBuf){
    	assert(isWriteB(lastBuf));
    	auto Mem_counterpart = lastBuf->write_other_ptr;	// Getting the Memory Write
    	assert(isWriteM(Mem_counterpart));
    	assert(!po.hasEdge(obsMem,Mem_counterpart));
		if(!po.hasEdge(Mem_counterpart,obsMem)){
			po.addEdge(Mem_counterpart,obsMem);		// Memory write before curr obs
		}
    }
  }
  // po.dump();
}

void ZClosure::preClose(const ZEvent *ev, const ZObs& obs){
  if (obs.thr == INT_MAX) {
  	// Handle initial-event observation separately
  	// No need to do for rule 1
  }else{
  	preClose(ev, po.basis.getEvent(obs));
  }
  // po.dump();
}
