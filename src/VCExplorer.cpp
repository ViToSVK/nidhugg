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

		auto nodesToMutate = getNodesToMutate();
    if (nodesToMutate.empty()) {
			// Fully executed trace
      ++executed_traces_full;
			current.reset();
			continue;
		}

		assert(current->unannot.empty());
		for (auto& nd : nodesToMutate)
			current->unannot.emplace(nd->getProcessID());
		
		// Get partial-order refinements that order extension writes
		// Each refinement will be a candidate for possible mutations
    std::list<PartialOrder> extendedPOs = extensionWritesOrderings();
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
/* GET NODES TO MUTATE         */
/* *************************** */

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

/* *************************** */
/* EXTENSION WRITES ORDERINGS  */
/* *************************** */

std::list<PartialOrder> VCExplorer::extensionWritesOrderings()
{
  assert(current.get());

	current->graph.initWorklist();

	if (current->annotation.empty()) {
		// No write is active yet, so nothing has to be ordered
		assert(current->graph.getExtensionFrom() == 0);
    return current->graph.dumpDoneWorklist();
	}

	// Order extension writes (same ml, different val):
	// active with active+nonactive, nonactive with active
	for (unsigned trace_idx = current->graph.getExtensionFrom();
			 trace_idx < current->trace.size(); ++trace_idx) {
		
    const VCEvent& ev = current->trace[trace_idx];
		if (isWrite(ev) && current->annotation.isActiveMl(ev)) {
      bool evIsActive = current->annotation.isActiveWrite(ev);
			current->graph.orderWrite(current->annotation, &ev, evIsActive);
		}
		
	}

	return current->graph.dumpDoneWorklist();
}

/* *************************** */
/* MUTATION ORDERINGS          */
/* *************************** */

std::list<PartialOrder> VCExplorer::mutationOrderings
(const PartialOrder& po, const Node *nd, int val,
 const VCAnnotation& annot, bool wasAlreadyActive)
{
  assert(current.get());

	// Close wrt annot including the new mutation 'nd sees val'
  bool initialclosure = current->graph.initWorklistAndClose(po, annot);


	if (!initialclosure) {
		// This mutation is impossible		
    return std::list<PartialOrder>();
	}
		
	
  if (wasAlreadyActive) {
		// This mutation does not activate any new writes
		// Our result is one unique mutated partial order in which
		// initialclosure made each read have only good tail writes
		// We get the above guarantee since
		// 'mutation activates nothing' -> 'everything active was
		// already fully ordered at the point of initialclosure'
    return current->graph.dumpDoneWorklist();
	}

	// Order all newly-active writes
	// TODO: RETAIN THESE FROM current->graph.getMutateValues ???
	// BUT IN THERE SOME EVENTS WE DISREGARDED SINCE HIDDEN, CAN WE EVEN GET AWAY WITHOUT ORDERING THESE HERE ???
	for (unsigned trace_idx = 0;
			 trace_idx < current->trace.size(); ++trace_idx) {
		
    const VCEvent& ev = current->trace[trace_idx];
		if (isWrite(ev) && ev.value == val && ev.ml == nd->getEvent()->ml) {
			// newly-active write
      assert(annot.isActiveWrite(ev)
						 && "We didn't properly update the annotation");
			current->graph.orderWrite(annot, &ev, true/*evIsActive*/);
		}
		
	}

	return current->graph.dumpDoneWorklist();
}

/* *************************** */
/* MUTATE READ                 */
/* *************************** */

void VCExplorer::mutateRead(const PartialOrder& po, const Node *nd)
{
  assert(isRead(nd->getEvent()));
	const VCEvent& ev = *(nd->getEvent());

	// Mutated 'threads with unannotated read'
	std::unordered_set<int> mutUnannot(current->unannot);
	assert(mutUnannot.count(nd->getProcessID()));
	mutUnannot.erase(nd->getProcessID());
	
	auto mutateValues = current->graph.getMutateValues(po, nd);

  /**/ llvm::errs() << "NODE[" << nd->getProcessID()
										<< "][" << nd->getEventID() << "] should see:";
	/**/ for (int val : mutateValues) llvm::errs() << " " << val;
	/**/ llvm::errs() << "\n\n";
	
	for (int val : mutateValues) {
		/**/ llvm::errs() << "MUTATING " << val << ", solutions:\n";

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

			/**/ current->graph.to_dot(mutPO, "");

			std::vector<VCEvent> mutTrace =
				extendTrace(current->graph.linearize(mutPO), mutUnannot);

      VCGraphVclock mutGraph(current->graph, /*base for graph*/
														 mutPO, /*base for 'original' po*/
														 mutTrace /*to extend the graph*/,
														 mutAnnotation /*to value-close extended 'original'*/);
			
			std::unique_ptr<VCTrace> mutatedVCTrace
				(new VCTrace(std::move(mutTrace),
										 mutAnnotation,
										 std::move(mutGraph)));

			worklist.push_back(std::move(mutatedVCTrace));
			assert(!mutatedVCTrace.get());
		}

	}
}

/* *************************** */
/* MUTATE LOCK                 */
/* *************************** */

void VCExplorer::mutateLock(const PartialOrder& po, const Node *nd)
{
  // TODO
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
