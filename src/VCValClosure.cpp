/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2018 Viktor Toman
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

#include "VCValClosure.h"

/* *************************** */
/* PREPARE ONE                 */
/* *************************** */

void VCValClosure::prepareOne
(const PartialOrder& po, const Node * readnd)
{
	assert(isRead(readnd->getEvent()));
	const SymAddrSize& ndml = readnd->getEvent()->ml;
	if (!wNonroot.count(ndml)) {
		// prepare wNonroot for a new ml
		wNonroot.emplace(ndml,
										 std::vector<const Node *>());
		auto itunord = graph.wNonrootUnord.find(ndml);
		assert(itunord != graph.wNonrootUnord.end());
		wNonroot[ndml].reserve(itunord->second.size());
		for (auto& nd : itunord->second)
			wNonroot[ndml].push_back(nd);
		std::sort(wNonroot[ndml].begin(), wNonroot[ndml].end(),
							VCGraphVclock::POcomp(graph, po));
	}

	// prepare wRem
	if (readnd->getProcessID() == graph.starRoot())
		wRem.emplace(readnd, wNonroot.at(ndml));
	else
		wRem.emplace(readnd, graph.wRoot.at(ndml));

	// prepare wBounds
	// Note: replace with binary searches?
	wBounds.emplace(readnd, std::pair<int, int>(0, wRem.at(readnd).size()-1));
	int& low = wBounds[readnd].first;
	int& high = wBounds[readnd].second;
	assert(low >= high || graph.hasEdge(wRem.at(readnd)[low],
																			wRem.at(readnd)[high], po));
	while (high >= 0) {
		if (graph.hasEdge(readnd, wRem.at(readnd)[high], po))
			--high; // this one is covered (since the read HB it)
		else
			break;
	}
	while (low < (int) wRem.at(readnd).size()) {
		if (low + 1 < (int) wRem.at(readnd).size() &&
				graph.hasEdge(wRem.at(readnd)[low + 1], readnd, po))
			++low; // this one is covered
		else
			break;
	}
	assert(low >= high || graph.hasEdge(wRem.at(readnd)[low],
																			wRem.at(readnd)[high], po));
	assert((high == -1 || !graph.hasEdge(readnd, wRem.at(readnd)[high], po))
				 && high + 1 <= (int) wRem.at(readnd).size() &&
				 (high + 1 == (int) wRem.at(readnd).size() ||
					graph.hasEdge(readnd, wRem.at(readnd)[high + 1], po))
				 && "badly computed high index into wRem (during prepare)");
	assert((low == (int) wRem.at(readnd).size()
					||
				 (
					graph.hasEdge(wRem.at(readnd)[low], readnd, po)
					&&
					(low + 1 == (int) wRem.at(readnd).size() ||
					 !graph.hasEdge(wRem.at(readnd)[low + 1], readnd, po))
				 )
					||
				 (
					!graph.hasEdge(wRem.at(readnd)[low], readnd, po) &&
					low == 0 // this also covers the case of read HB first write
				 )) && "badly computed low index into wRem (during prepare)");

	// prepare wLoc
	// Note: for a non-root read, 'local' means
	// from any non-root process
	const std::vector<const Node *>& localwrites =
		(readnd->getProcessID() == graph.starRoot())
		? graph.wRoot.at(ndml)
		: wNonroot.at(ndml);

	int idx = localwrites.size() - 1;
	while (idx >= 0) {
		if (!graph.hasEdge(localwrites[idx], readnd, po)) {
			assert(graph.hasEdge(readnd, localwrites[idx], po)
						 && "These should be fully ordered");
			--idx;
		} else
			break;
	}

  const Node * localcandidate =
		(idx >= 0) ? localwrites[idx] : graph.initial_node;

  assert(graph.hasEdge(localcandidate, readnd, po));
	
	// Check if the local write is covered
	// It is sufficient to check if it is covered
	// only by the 'low' remote write since all
	// remote writes after 'low' do not HB readnd
	// (otherwise 'low' itself would be covered)
	if (low < (int) wRem.at(readnd).size() &&
			graph.hasEdge(wRem.at(readnd)[low], readnd, po) &&
			graph.hasEdge(localcandidate, wRem.at(readnd)[low], po))
		wLoc.emplace(readnd, nullptr); // local write is covered
	else {
		wLoc.emplace(readnd, localcandidate); // local write is not covered
		// now check if the local write covers
		// the 'low' remote write
		if (low < (int) wRem.at(readnd).size() &&
				graph.hasEdge(wRem.at(readnd)[low], localcandidate, po)) {
			// it does
      assert(graph.hasEdge(wRem.at(readnd)[low], readnd, po));
			++low;
			assert(low == (int) wRem.at(readnd).size() ||
						 !graph.hasEdge(wRem.at(readnd)[low], readnd, po));
		}
	}

}

/* *************************** */
/* PREPARE                     */
/* *************************** */

void VCValClosure::prepare
(const PartialOrder& po, const Node * newread)
{
  if (newread) {
    // performing a new annotation
		// caches for all old annotations are set
		#ifndef NDEBUG
		for (auto& fun : valFunction) {
      assert(wRem.count(fun.first) && wBounds.count(fun.first)
						 && wLoc.count(fun.first));
		}
		#endif

		prepareOne(po, newread);
	} else {
		// either checking an extension against
		// previous annotation, or closing after
		// a lock annotation
		if (!wRem.empty()) {
			#ifndef NDEBUG
			for (auto& fun : valFunction) {
				assert(wRem.count(fun.first) && wBounds.count(fun.first)
							 && wLoc.count(fun.first));
			}
			#endif
		} else {
			for (auto& fun : valFunction)
				prepareOne(po, fun.first);
		}

	}
}

/* *************************** */
/* UPDATE VISIBLE CACHE        */
/* *************************** */

void VCValClosure::updateVisibleCache
(const PartialOrder& po, const Node * readnd)
{
	assert(wRem.count(readnd) && wBounds.count(readnd)
				 && wLoc.count(readnd));

	// Update remote
	int& low = wBounds[readnd].first;
	int& high = wBounds[readnd].second;
	assert(low >= high || graph.hasEdge(wRem.at(readnd)[low],
																			wRem.at(readnd)[high], po));
	while (high >= 0) {
		if (graph.hasEdge(readnd, wRem.at(readnd)[high], po))
			--high; // this one is covered (since the read HB it)
		else
			break;
	}
	while (low < (int) wRem.at(readnd).size()) {
		if (low + 1 < (int) wRem.at(readnd).size() &&
				graph.hasEdge(wRem.at(readnd)[low + 1], readnd, po))
			++low; // this one is covered
		else
			break;
	}
	assert(low >= high || graph.hasEdge(wRem.at(readnd)[low],
																			wRem.at(readnd)[high], po));
	assert((high == -1 || !graph.hasEdge(readnd, wRem.at(readnd)[high], po))
				 && high + 1 <= (int) wRem.at(readnd).size() &&
				 (high + 1 == (int) wRem.at(readnd).size() ||
					graph.hasEdge(readnd, wRem.at(readnd)[high + 1], po))
				 && "badly computed high index into wRem (during update)");
	assert((low == (int) wRem.at(readnd).size()
					||
				 (
					graph.hasEdge(wRem.at(readnd)[low], readnd, po)
					&&
					(low + 1 == (int) wRem.at(readnd).size() ||
					 !graph.hasEdge(wRem.at(readnd)[low + 1], readnd, po))
				 )
					||
				 (
					!graph.hasEdge(wRem.at(readnd)[low], readnd, po) &&
					low == 0 // this also covers the case of read HB first write
				 )) && "badly computed low index into wRem (during update)");

	// Update local
	if (wLoc[readnd]) {
		// Check if the local write is covered
		// It is sufficient to check if it is covered
		// only by the 'low' remote write since all
		// remote writes after 'low' do not HB readnd
		// (otherwise 'low' itself would be covered)
		if (low < (int) wRem.at(readnd).size() &&
				graph.hasEdge(wRem.at(readnd)[low], readnd, po) &&
				graph.hasEdge(wLoc[readnd], wRem.at(readnd)[low], po))
			wLoc[readnd] = nullptr; // local write is covered
		else {
			// local write is not covered
			// now check if the local write covers
			// the 'low' remote write
			if (low < (int) wRem.at(readnd).size() &&
					graph.hasEdge(wRem.at(readnd)[low], wLoc[readnd], po)) {
				// it does
				assert(graph.hasEdge(wRem.at(readnd)[low], readnd, po));
				++low;
				assert(low == (int) wRem.at(readnd).size() ||
							 !graph.hasEdge(wRem.at(readnd)[low], readnd, po));
			}
		}
	}
	
}

/* *************************** */
/* RULE 1                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> VCValClosure::ruleOne
(const PartialOrder& po, const Node * readnd,
 const std::pair<int, VCAnnotation::Loc>& val)
{
	assert(isRead(readnd->getEvent()));
	// Rule1: there exists a write such that
	// (i) GOOD (ii) HEAD (iii) HB read
	//
	// Note: the local write is always
	// either *HEAD* or *covered*

  if (readnd->getProcessID() != graph.starRoot()) {
		// Non-root readnd with ANY location
    assert(val.second == VCAnnotation::Loc::ANY);

		return {false, false};
	}

	assert(readnd->getProcessID() == graph.starRoot());

	if (val.second == VCAnnotation::Loc::LOCAL) {
    // Root readnd with LOCAL location

		return {false, false};
	}

	assert(val.second == VCAnnotation::Loc::REMOTE);

	// Root readnd with REMOTE location

	return {false, false};
}

/* *************************** */
/* RULE 2                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> VCValClosure::ruleTwo
(const PartialOrder& po, const Node * readnd,
 const std::pair<int, VCAnnotation::Loc>& val)
{
	assert(isRead(readnd->getEvent()));
	// Rule2: there exists a write such that
	// (i) GOOD (ii) TAIL
	//
	// Note: if there is a local and a remote
	// GOOD write, make TAIL the one which
	// forces to add a weaker edge

  if (readnd->getProcessID() != graph.starRoot()) {
		// Non-root readnd with ANY location
    assert(val.second == VCAnnotation::Loc::ANY);

		return {false, false};
	}

	assert(readnd->getProcessID() == graph.starRoot());

	if (val.second == VCAnnotation::Loc::LOCAL) {
    // Root readnd with LOCAL location

		return {false, false};
	}

	assert(val.second == VCAnnotation::Loc::REMOTE);

	// Root readnd with REMOTE location

	return {false, false};
}

/* *************************** */
/* RULE 3                      */
/* *************************** */

// first) true iff impossible
// second) true iff something changed
std::pair<bool, bool> VCValClosure::ruleThree
(const PartialOrder& po, const Node * readnd,
 const std::pair<int, VCAnnotation::Loc>& val)
{
	assert(isRead(readnd->getEvent()));
	// Rule3: for every write such that
	// (i) BAD (ii) HEAD (iii) HB read
	// it HB some visible GOOD one
	//
	// Note: after Rule1 there is at most
	// one write satisfying (i,ii,iii), and
	// if it doesn't HB some visible GOOD one,
	// we need it to HB the 'other' tail
	// (which is GOOD after Rule2)

  if (readnd->getProcessID() != graph.starRoot()) {
		// Non-root readnd with ANY location
    assert(val.second == VCAnnotation::Loc::ANY);

		return {false, false};
	}

	assert(readnd->getProcessID() == graph.starRoot());

	if (val.second == VCAnnotation::Loc::LOCAL) {
    // Root readnd with LOCAL location

		return {false, false};
	}

	assert(val.second == VCAnnotation::Loc::REMOTE);

	// Root readnd with REMOTE location

	return {false, false};
}

/* *************************** */
/* VAL CLOSE                   */
/* *************************** */

void VCValClosure::valClose
(const PartialOrder& po, const Node * newread,
 const std::pair<int, VCAnnotation::Loc> * newval)
{
  prepare(po, newread);
	
  bool change = true;
	while (change) {
    change = false;
		for (auto& read_val : valFunction) {
      // Update the visible writes cache
			updateVisibleCache(po, read_val.first);
      // Rule1
      auto res = ruleOne(po, read_val.first, read_val.second);
			if (res.first) { closed = false; return; }
			if (res.second) change = true;
			// Rule2
      res = ruleTwo(po, read_val.first, read_val.second);
			if (res.first) { closed = false; return; }
			if (res.second) change = true;
      //Rule3
      res = ruleThree(po, read_val.first, read_val.second);
			if (res.first) { closed = false; return; }
			if (res.second) change = true;
		}
		if (newread) {
      // Update the visible writes cache
			updateVisibleCache(po, newread);
      // Rule1
      auto res = ruleOne(po, newread, *newval);
			if (res.first) { closed = false; return; }
			if (res.second) change = true;
			// Rule2
      res = ruleTwo(po, newread, *newval);
			if (res.first) { closed = false; return; }
			if (res.second) change = true;
      //Rule3
      res = ruleThree(po, newread, *newval);
			if (res.first) { closed = false; return; }
			if (res.second) change = true;
		}
	}

	assert(!change);
	closed = true;
}

/* *************************** */
/* VAL CLOSE LOCK              */
/* *************************** */

void VCValClosure::valCloseLock(const PartialOrder& po,
																const Node * locknode,
																const Node * lastunlocknode)
{
  assert(isLock(locknode->getEvent()) &&
				 isUnlock(lastunlocknode->getEvent()) &&
				 locknode->getEvent()->ml == lastunlocknode->getEvent()->ml);

	if (graph.hasEdge(locknode, lastunlocknode, po)) {
    closed = false;
		return;
	} else if (graph.hasEdge(lastunlocknode, locknode, po)) {
    closed = true;
		return;
	}

	graph.addEdge(lastunlocknode, locknode, po);
	
	valClose(po, nullptr, nullptr);
}
