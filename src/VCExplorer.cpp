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
	llvm::dbgs() << "\n!!!DONE!!!\n";	
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
		/*
    llvm::errs() << "********* TRACE *********\n";                               ///////////////////////
		current->annotation.dump();
		current->graph.to_dot("");
		// if (current->graph.getExtensionFrom() > 0) return;
		*/
		auto nodesToMutate = current->graph.getNodesToMutate();
		for (auto it = nodesToMutate.begin(); it != nodesToMutate.end(); ) {
			const Node * nd = *it;
			assert(isRead(nd->getEvent()) || isLock(nd->getEvent()));
			if ((isRead(nd->getEvent()) && current->annotation.defines(nd)) ||
					(isLock(nd->getEvent()) && current->annotation.isLastLock(nd)))
        it = nodesToMutate.erase(it);
      else
        ++it;
    }
    if (nodesToMutate.empty()) {
			// Fully executed trace
      ++executed_traces_full;
			//llvm::errs() << "********* FULL TRACE *********\n";                        ///////////////////////
			//current->annotation.dump();
			//current->graph.to_dot("");
			current.reset();
			continue;
		}

		assert(current->unannot.empty());
		for (auto& nd : nodesToMutate)
			current->unannot.insert(nd->getEvent()->iid.get_pid());
		
		// Get partial-order refinements that order extension events
		// Each refinement will be a candidate for possible mutations
    std::list<PartialOrder> extendedPOs = extensionEventsOrderings();
		assert(!extendedPOs.empty() && "current->trace is one witness");
		
    while (!extendedPOs.empty()) {
      auto po = std::move(extendedPOs.front());
			assert(!extendedPOs.front().first.get() &&
						 !extendedPOs.front().second.get());
			extendedPOs.pop_front();

			/*
			llvm::errs() << "********* EXTENSION *********\n";                           ///////////////////////
			current->graph.to_dot(po, "");
			*/
			
			auto withoutMutation = VCValClosure(current->graph, current->annotation);
			withoutMutation.valClose(po, nullptr, nullptr);
			if (!withoutMutation.closed) {
        // This ordering of extension events
				// makes the original annotation unrealizable
				// therefore no need to try any mutations
				llvm::errs() << "Invalidates annotation we have so far\n";                 // ///////////////////////
				po.first.reset();
				po.second.reset();
				continue;
			}

			// Try all possible mutations
			for (auto nd : nodesToMutate) {
        if (isRead(nd->getEvent()))
					mutateRead(po, withoutMutation, nd);
				else {
          assert(isLock(nd->getEvent()));
					mutateLock(po, withoutMutation, nd);
				}
			}

			po.first.reset();
			po.second.reset();
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

void VCExplorer::mutateRead(const PartialOrder& po, const VCValClosure& withoutMutation, const Node *nd)
{
  assert(isRead(nd->getEvent()));
	
	std::unordered_set<int> mutatedUnannot(current->unannot);
	assert(mutatedUnannot.count(nd->getEvent()->iid.get_pid()));
	mutatedUnannot.erase(nd->getEvent()->iid.get_pid());
	
	auto mutateValues =
		current->graph.getMutateValues(po, nd);

	/*
  llvm::errs() << "NODE[" << nd->getProcessID()
							 << "][" << nd->getEventID() << "] should see:";                      ///////////////////////
	for (auto& val_pos : mutateValues) {
    llvm::errs() << " " << val_pos.first << "-";
		char loc = (val_pos.second == VCAnnotation::Loc::LOCAL)?'L':
			((val_pos.second == VCAnnotation::Loc::REMOTE)?'R':'A');
		llvm::errs() << loc;
	}
	llvm::errs() << "\n\n";
	*/

	for (auto& val_pos : mutateValues) {
		/*
		llvm::errs() << val_pos.first << "-";
		char loc = (val_pos.second == VCAnnotation::Loc::LOCAL)?'L':                     ///////////////////////
			((val_pos.second == VCAnnotation::Loc::REMOTE)?'R':'A');
		llvm::errs() << loc << "... ";
		*/
		auto mutatedPo = PartialOrder(std::unique_ptr<ThreadPairsVclocks>
																	(new ThreadPairsVclocks(*(po.first))),
																	std::unique_ptr<ThreadPairsVclocks>
																	(new ThreadPairsVclocks(*(po.second))));
		
    auto withMutation = VCValClosure(withoutMutation);
		// withMutation.valClose(mutatedPo, nd, &val_pos);                  TODO
		
		if (!withMutation.closed) {
			// The mutation on 'mutatedPo' failed
		  // llvm::errs() << "FAILED\n";                                                         ///////////////////////
      // current->graph.to_dot(mutatedPo, "");			
			mutatedPo.first.reset();
			mutatedPo.second.reset();
			continue;
		}

		// The mutation on 'mutatedPo' succeeded
		// llvm::errs() << "SUCCEEDED\n";                                                    ///////////////////////
    // current->graph.to_dot(mutatedPo, "");
		
		VCAnnotation mutatedAnnotation(current->annotation);
		// mutatedAnnotation.add(*(nd->getEvent()), val_pos);                  TODO

		std::vector<VCEvent> mutatedTrace =
			extendTrace(current->graph.linearize(mutatedPo), mutatedUnannot);
		assert(traceRespectsAnnotation(mutatedTrace, mutatedAnnotation));

		VCGraphVclock mutatedGraph(current->graph,       // base for graph
															 std::move(mutatedPo), // base for 'original' po
														   mutatedTrace);        // to extend the graph
    assert(!mutatedPo.first.get() && !mutatedPo.second.get());

		std::unique_ptr<VCTrace> mutatedVCTrace
			(new VCTrace(std::move(mutatedTrace),
									 std::move(mutatedAnnotation),
									 std::move(mutatedGraph)));
		assert(mutatedTrace.empty() && mutatedAnnotation.empty() && mutatedGraph.empty());

		worklist.push_back(std::move(mutatedVCTrace));
		assert(!mutatedVCTrace.get());
	
	}
	
}

/* *************************** */
/* MUTATE LOCK                 */
/* *************************** */

void VCExplorer::mutateLock(const PartialOrder& po, const VCValClosure& withoutMutation, const Node *nd)
{
  assert(isLock(nd->getEvent()));
	
  auto lastLock = current->annotation.getLastLock(nd);

  if (!lastLock.first) {
    // This lock hasn't been touched before
		// Trivially realizable

		std::unordered_set<int> mutatedUnannot(current->unannot);
		assert(mutatedUnannot.count(nd->getEvent()->iid.get_pid()));
		mutatedUnannot.erase(nd->getEvent()->iid.get_pid());
		
		auto mutatedPo = PartialOrder(std::unique_ptr<ThreadPairsVclocks>
																	(new ThreadPairsVclocks(*(po.first))),
																	std::unique_ptr<ThreadPairsVclocks>
																	(new ThreadPairsVclocks(*(po.second))));
		
		VCAnnotation mutatedAnnotation(current->annotation);
		mutatedAnnotation.setLastLock(nd);

		std::vector<VCEvent> mutatedTrace =
			extendTrace(current->graph.linearize(mutatedPo), mutatedUnannot);
		assert(traceRespectsAnnotation(mutatedTrace, mutatedAnnotation));

		VCGraphVclock mutatedGraph(current->graph,       // base for graph
															 std::move(mutatedPo), // base for 'original' po
														   mutatedTrace);        // to extend the graph
    assert(!mutatedPo.first.get() && !mutatedPo.second.get());

		std::unique_ptr<VCTrace> mutatedVCTrace
			(new VCTrace(std::move(mutatedTrace),
									 std::move(mutatedAnnotation),
									 std::move(mutatedGraph)));
		assert(mutatedTrace.empty() && mutatedAnnotation.empty() && mutatedGraph.empty());

		worklist.push_back(std::move(mutatedVCTrace));
		assert(!mutatedVCTrace.get());
		return;
	}

	assert(lastLock.first);
	auto lastunlockit = current->graph.nodes_iterator(lastLock.second);
	#ifndef NDEBUG
  const Node * lastlocknd = *lastunlockit;
  #endif
	const Node * lastunlocknd = nullptr;
	
	while (!lastunlockit.atProcessEnd()) {
    ++lastunlockit;
		const VCEvent& cand = *((*lastunlockit)->getEvent());
		if (isUnlock(cand) && cand.ml == nd->getEvent()->ml) {
      assert(cand.ml == lastlocknd->getEvent()->ml);
			lastunlocknd = (*lastunlockit);
			break;
		}
	}

	if (!lastunlocknd) {
    // This lock is currently locked
		// Trivially unrealizable
		return;
	}

	assert(lastunlocknd);

	auto mutatedPo = PartialOrder(std::unique_ptr<ThreadPairsVclocks>
																(new ThreadPairsVclocks(*(po.first))),
																std::unique_ptr<ThreadPairsVclocks>
																(new ThreadPairsVclocks(*(po.second))));

	auto withMutation = VCValClosure(withoutMutation);
	withMutation.valCloseLock(mutatedPo, nd, lastunlocknd);

	if (!withMutation.closed) {
		// The mutation on 'mutatedPo' failed
		mutatedPo.first.reset();
		mutatedPo.second.reset();
		return;
	}

	// The mutation on 'mutatedPo' succeeded

	std::unordered_set<int> mutatedUnannot(current->unannot);
	assert(mutatedUnannot.count(nd->getEvent()->iid.get_pid()));
	mutatedUnannot.erase(nd->getEvent()->iid.get_pid());
	
	VCAnnotation mutatedAnnotation(current->annotation);
	mutatedAnnotation.setLastLock(nd);

	std::vector<VCEvent> mutatedTrace =
		extendTrace(current->graph.linearize(mutatedPo), mutatedUnannot);
	assert(traceRespectsAnnotation(mutatedTrace, mutatedAnnotation));

	VCGraphVclock mutatedGraph(current->graph,       // base for graph
														 std::move(mutatedPo), // base for 'original' po
														 mutatedTrace);        // to extend the graph
	assert(!mutatedPo.first.get() && !mutatedPo.second.get());

	std::unique_ptr<VCTrace> mutatedVCTrace
		(new VCTrace(std::move(mutatedTrace),
								 std::move(mutatedAnnotation),
								 std::move(mutatedGraph)));	
	assert(mutatedTrace.empty() && mutatedAnnotation.empty() && mutatedGraph.empty());

	worklist.push_back(std::move(mutatedVCTrace));
	assert(!mutatedVCTrace.get());

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

/* *************************** */
/* TRACE RESPECTS ANNOTATION   */
/* *************************** */

bool VCExplorer::traceRespectsAnnotation(const std::vector<VCEvent>& trace,
																				 const VCAnnotation& annotation) const {
	// TODO: check also if observes one of the good writes
  for (unsigned i=0; i < trace.size(); ++i) {
		const VCEvent& ev = trace[i];
		if (isRead(ev) && annotation.defines(ev.pid, ev.event_order)) {
      const auto& ann = annotation.getAnn(ev.pid, ev.event_order);
			if (ann.value != ev.value)
				return false;
			for (int j=i-1; j >= -1; --j) {
        if (j == -1) {
					
          if (ann.value != 0 || ev.value != 0)
						return false;
					if ((int) current->graph.starRoot() == ev.iid.get_pid() / 2 &&
							ann.loc != VCAnnotation::Loc::LOCAL)
						return false;
					if ((int) current->graph.starRoot() != ev.iid.get_pid() / 2 &&
							ann.loc != VCAnnotation::Loc::ANY)
						return false;
					break;
					
				} else {
					
          const VCEvent& wrev = trace[j];
          if (isWrite(wrev) && wrev.ml == ev.ml) {
						assert(wrev.value == ev.value);
						if (ann.value != wrev.value)
							return false;
						if ((int) current->graph.starRoot() == ev.iid.get_pid() / 2 &&
								ev.iid.get_pid() == wrev.iid.get_pid() &&
								ann.loc != VCAnnotation::Loc::LOCAL)
							return false;
						if ((int) current->graph.starRoot() == ev.iid.get_pid() / 2 &&
								ev.iid.get_pid() != wrev.iid.get_pid() &&
								ann.loc != VCAnnotation::Loc::REMOTE)
							return false;						
						if ((int) current->graph.starRoot() != ev.iid.get_pid() / 2 &&
								ann.loc != VCAnnotation::Loc::ANY)
							return false;
						break;
					}
					
				}
			}
		}
	}

	return true;
}
