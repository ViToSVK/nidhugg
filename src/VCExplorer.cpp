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

#include <iostream>

#include "VCExplorer.h"
#include "VCHelpers.h"

void VCExplorer::print_stats()
{
  llvm::dbgs() << "Executed traces: " << executed_traces << "\n";
  //llvm::dbgs() << "Total number of executed instructions: " << total_instr_executed << "\n";

	// Change to false to test if assertions are on
	// To disable assertions (i.e. build as Release),
	// in src/Makefile add in CXXFLAGS this: -DNDEBUG
	assert(true && "TEST");
}

void VCExplorer::explore()
{
  while (!worklist.empty()) {

		// Get a VCTrace
		assert(!current.get());
		current = std::move(worklist.front());
		assert(!worklist.front().get());
		worklist.pop_front();

		// For a given memory location, have a map
		// read -> values written into that location
		std::unordered_map<SymAddrSize,
											 std::unordered_map<const Node *, std::set<int>>>
			reads_to_annotate;

		// For a given memory location, have a map
		// value -> writes that wrote it into that location
		std::unordered_map<SymAddrSize,
											 std::unordered_map<int, std::set<const Node *>>>
			active_writes;


		// Find the reads to annotate
		for (unsigned proc_idx = 0;
				 proc_idx < current->partialOrder.size();
				 ++proc_idx) {

			auto it = current->partialOrder.nodes_process_end(proc_idx);
      const Node *nd = *it;
			const VCEvent *ev = nd->getEvent();
			
			if (ev->may_conflict && isRead(ev) &&
					!current->annotation.defines(*ev) ) {
        auto it = reads_to_annotate.find(ev->ml);
				if (it == reads_to_annotate.end()) {
					reads_to_annotate.emplace_hint(it, ev->ml,
																				 std::unordered_map<const Node *, std::set<int>>());
				}
				reads_to_annotate[ev->ml].emplace(nd, std::set<int>());
				reads_to_annotate[ev->ml][nd].insert(0); // initial node XXX FOR SPAWNED THREADS ALSO, RIGHT?
			}
		}

		// Find the active writes
		

		
			

		// Delete managed VCTrace
		current.reset();
	}
}

std::vector<VCEvent> VCExplorer::extendTrace
(std::vector<VCEvent>&& tr, const std::unordered_set<int>& unannot) {
	VCTraceBuilder TB(originalTB.config, originalTB.M, std::move(tr), unannot);
	auto res = TB.extendGivenTrace();

	if (TB.has_error())
		originalTB.error_trace = TB.get_trace();

	++executed_traces;

	return res;
}
