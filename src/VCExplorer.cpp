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
	llvm::dbgs() << "Fully executed traces:            " << executed_traces_full << "\n";
	llvm::dbgs() << "Fully+partially executed traces:  " << executed_traces << "\n";
	llvm::dbgs() << "\n";

	// Change to false to test if assertions are on
	// To disable assertions (i.e. build as Release),
	// in src/Makefile add in CXXFLAGS this: -DNDEBUG
	assert(true && "RUN ON RELEASE");
}

/* *************************** */
/* EXPLORE                     */
/* *************************** */

void VCExplorer::explore()
{
  while (!worklist.empty()) {

		// Get a VCTrace
		assert(!current.get());
		current = std::move(worklist.front());
		assert(!worklist.front().get());
		worklist.pop_front();

		auto nodesToMutate = current->graph.getNodesToMutate();
    if (nodesToMutate.empty()) {
			// Fully executed trace
      ++executed_traces_full;
			current.reset();
			continue;
		}

		assert(current->unannot.empty());
		for (auto& nd : nodesToMutate)
			current->unannot.insert(nd->getProcessID());
		
		// Get partial-order refinements that order extension events
		// Each refinement will be a candidate for possible mutations
    std::list<PartialOrder> extendedPOs = extensionEventsOrderings();
		assert(!extendedPOs.empty() && "current->trace is one witness");

		/**/ if (current->graph.getExtensionFrom() > 20) return;
		/**/ llvm::errs() << "********* TRACE *********\n";
		
    while (!extendedPOs.empty()) {
      auto po = std::move(extendedPOs.front());
			assert(!extendedPOs.front().first.get());
			assert(!extendedPOs.front().second.get());
			extendedPOs.pop_front();

			/**/ llvm::errs() << "********* EXTENSION *********\n";
			/**/ current->graph.to_dot(po, "");

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

/* *************************** */
/* EXTENSION EVENTS ORDERINGS  */
/* *************************** */

std::list<PartialOrder> VCExplorer::extensionEventsOrderings()
{
  assert(current.get());

	current->graph.initWorklist();
	
	// Order extension reads + writes of non-star processes
	// with conflicting writes of non-star processes
	for (unsigned trace_idx = current->graph.getExtensionFrom();
			 trace_idx < current->trace.size(); ++trace_idx) {
		
    const VCEvent& ev = current->trace[trace_idx];
		const Node *nd = current->graph.getNode(ev);
		if ((isWrite(ev) || isRead(ev)) && nd->getProcessID() != current->graph.starRoot()) {
			current->graph.orderEventMaz(&ev);
		}
		
	}

	return current->graph.dumpDoneWorklist();
}

/* *************************** */
/* MUTATE READ                 */
/* *************************** */

void VCExplorer::mutateRead(const PartialOrder& po, const Node *nd)
{
	const VCEvent& ev = *(nd->getEvent());
  assert(isRead(ev));
	bool isFromStarRoot =
		(nd->getProcessID() == current->graph.starRoot());
	
	// Mutated 'threads with unannotated read'
	std::unordered_set<int> mutUnannot(current->unannot);
	assert(mutUnannot.count(nd->getProcessID()));
	mutUnannot.erase(nd->getProcessID());
	
	auto mutateValues =
		current->graph.getMutateValues(po, nd); // , isFromStarRoot

  /**/ llvm::errs() << "NODE[" << nd->getProcessID()
										<< "][" << nd->getEventID() << "] should see:";
	/**/ for (int val : mutateValues) llvm::errs() << " " << val;
	/**/ llvm::errs() << "\n\n";

	/*
	for (int val : mutateValues) {
		llvm::errs() << "MUTATING " << val << ", solutions:\n"; ///////

		// Mutated annotation
		VCAnnotation mutAnnotation(current->annotation);
		mutAnnotation.add(ev, val);
		bool wasAlreadyActive = mutAnnotation.isActiveMlVal(ev, val);
		if (!wasAlreadyActive)
			mutAnnotation.addActiveValue(ev, val);

		// Get partial-order refinements that successfully mutate 'po'
		// such that read event-node 'nd' observes value 'val'
		std::list<PartialOrder> mutatedPOs =
			mutationOrderings(po, nd, val, mutAnnotation, wasAlreadyActive);

		// For each mutated partial order:
	  // 1) linearize it into a trace 2) extend this linearization,
		// 3) create a new corresponding graph 4) add to worklist
    while (!mutatedPOs.empty()) {
      auto mutPO = std::move(mutatedPOs.front());
			assert(!mutatedPOs.front().first.get());
			assert(!mutatedPOs.front().second.get());
			mutatedPOs.pop_front();

			current->graph.to_dot(mutPO, ""); ///////

			std::vector<VCEvent> mutTrace =
				extendTrace(current->graph.linearize(mutPO), mutUnannot);

      VCGraphVclock mutGraph(current->graph, // base for graph
														 mutPO, // base for 'original' po
														 mutTrace, // to extend the graph
														 mutAnnotation); // to value-close extended 'original'
      
			std::unique_ptr<VCTrace> mutatedVCTrace
				(new VCTrace(std::move(mutTrace),
										 mutAnnotation,
										 std::move(mutGraph)));
			assert(mutTrace.empty() && mutGraph.empty());

			worklist.push_back(std::move(mutatedVCTrace));
			assert(!mutatedVCTrace.get());
		}

	}
  */
}

/* *************************** */
/* MUTATE LOCK                 */
/* *************************** */

void VCExplorer::mutateLock(const PartialOrder& po, const Node *nd)
{
	const VCEvent& ev = *(nd->getEvent());
  assert(isLock(ev));

	// Mutated 'threads with unannotated read'
	std::unordered_set<int> mutUnannot(current->unannot);
	assert(mutUnannot.count(nd->getProcessID()));
	mutUnannot.erase(nd->getProcessID());

	// TODO CONTINUE
	
}

/* *************************** */
/* EXTEND TRACE                */
/* *************************** */

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
