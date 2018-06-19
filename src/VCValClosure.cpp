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

/* *************************** */
/* VALUE CLOSURE               */
/* *************************** */

#include "VCValClosure.h"

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
			--high; // since the read HB 
		else
			break;
	}
	while (low < (int) wRem.at(readnd).size()) {
		if (low + 1 < (int) wRem.at(readnd).size() &&
				graph.hasEdge(wRem.at(readnd)[low + 1], readnd, po))
			++low; // since this one is covered
		else
			break;
	}
	assert(low >= high || graph.hasEdge(wRem.at(readnd)[low],
																			wRem.at(readnd)[high], po));
	assert((high == -1 || high == (int) wRem.at(readnd).size()-1 ||
					graph.hasEdge(readnd, wRem.at(readnd)[high], po))
				 && (high - 1 < 0 ||
						 !graph.hasEdge(readnd, wRem.at(readnd)[high - 1], po))
				 && "badly computed high index into wRem");
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
					low == 0
				 )) && "badly computed low index into wRem");

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

	if (idx >= 0)
		wLoc.emplace(readnd, localwrites[idx]);
	else
		wLoc.emplace(readnd, graph.initial_node);
}

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
    for (auto& fun : valFunction)
      prepareOne(po, fun.first);
	}
}

void VCValClosure::valClose
(const PartialOrder& po, const Node * newread,
 const std::pair<int, VCAnnotation::Loc> * newval)
{
  bool change = true;
	while (change) {
    change = false;
		for (auto& fun : valFunction) {
      // call rule 1

			// call rule 2

			// call rule 3
		}

		if (newread) {
      // call rule 1

			// call rule 2

			// call rule 3
		}
		
	}
}
