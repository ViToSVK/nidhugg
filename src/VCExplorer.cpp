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
  llvm::dbgs() << "Total number of executed instructions: " << total_instr_executed << "\n";

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

		// Get partial-order refinements that order extension writes
		// Each refinement will be a candidate for possible mutations
    std::list<PartialOrder> mutationPOs = extensionWritesOrderings();
		
    for (auto it = mutationPOs.begin();
				 it != mutationPOs.end(); ++it)
		  current->graph.to_dot(*it, "");

		/**/ if (current->graph.getExtensionFrom() > 0) return;

		auto nodesToMutate = getNodesToMutate();
		auto unannot = std::unordered_set<int>();
		for (auto& nd : nodesToMutate)
			unannot.emplace(nd->getProcessID());
		
    while (!mutationPOs.empty()) {
      auto po = std::move(mutationPOs.front());
			assert(!mutationPOs.front().first.get());
			assert(!mutationPOs.front().second.get());
			mutationPOs.pop_front();

			for (auto nd : nodesToMutate) {
        if (isRead(nd->getEvent()))
					mutateRead(po, nd);
				else {
          assert(isLock(nd->getEvent()));
					mutateLock(po, nd);
				}
			}
		}
		
		// Delete managed VCTrace
		current.reset();
	}
}

std::list<PartialOrder> VCExplorer::extensionWritesOrderings()
{
  assert(current.get());

	current->graph.initWorklist();

	if (current->annotation.empty()) {
		// No write is active yet, so nothing has to be ordered
		assert(current->graph.getExtensionFrom() == 0);
    return current->graph.dumpMutationCandidates();
	}

	for (unsigned trace_idx = current->graph.getExtensionFrom();
			 trace_idx < current->trace.size(); ++trace_idx) {
		
    const VCEvent& ev = current->trace[trace_idx];
		if (isWrite(ev) && current->annotation.isActiveMl(ev)) {
      bool evIsActive = current->annotation.isActiveWrite(ev);
			current->graph.orderWrite(current->annotation, &ev, evIsActive);
		}
		
	}

	return current->graph.dumpMutationCandidates();
}

std::vector<VCEvent> VCExplorer::extendTrace(std::vector<VCEvent>&& tr,
																						 const std::unordered_set<int>& unannot)
{
	VCTraceBuilder TB(originalTB.config, originalTB.M, std::move(tr), unannot);
	auto res = TB.extendGivenTrace();

	if (TB.has_error())
		originalTB.error_trace = TB.get_trace();

	++executed_traces;

	return res;
}

std::unordered_set<const Node *> VCExplorer::getNodesToMutate()
{
  auto candidates = current->graph.getLastNodes();
	
	auto result = std::unordered_set<const Node *>();

	for (auto& nd : candidates)
		if ((isRead(nd->getEvent()) || isLock(nd->getEvent()))
				&& !current->annotation.defines(*(nd->getEvent())))
			result.emplace(nd);

	return result;
}

void VCExplorer::mutateRead(const PartialOrder& po, const Node *nd)
{
  assert(isRead(nd->getEvent()));

	auto mutateValues = current->graph.getMutateValues(po, nd);

	for (int val : mutateValues) {
    // mutate read-value
	}
}

void VCExplorer::mutateLock(const PartialOrder& po, const Node *nd)
{
  // TODO
}
