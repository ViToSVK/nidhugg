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
#include <iomanip>

#include "VCExplorer.h"
#include "VCHelpers.h"

void VCExplorer::print_stats()
{
	std::setprecision(4);
	std::cout << "\n";	
	std::cout << "Fully executed traces:            " << executed_traces_full << "\n";
	std::cout << "Fully+partially executed traces:  " << executed_traces << "\n";
	std::cout << "Total time spent on copying:      " << time_graphcopy << "\n";
	std::cout << "Total time spent on replaying:    " << time_replaying << "\n";
	std::cout << "Total time spent on mazurkiewicz: " << time_maz << "\n";	
	std::cout << "Total time spent on closure:      " << time_closure << "\n";
	std::cout << "Closure after ordering failed:    " << cl_ordering_failed << "\n";
	std::cout << "Closure after ordering succeeded: " << cl_ordering_succeeded << "\n";
	std::cout << "Closure after mutation failed:    " << cl_mutation_failed << "\n";
	std::cout << "Closure after mutation succeeded: " << cl_mutation_succeeded << "\n";
	std::cout << "\n";

	// Change to false to test if assertions are on
	// To disable assertions (i.e. build as Release),
	// in src/Makefile add in CXXFLAGS this: -DNDEBUG
	assert(true && "RUN ON RELEASE");
}

/* *************************** */
/* EXPLORE                     */
/* *************************** */

bool VCExplorer::explore()
{
  while (!worklist.empty()) {
		// Get a VCTrace
		assert(!current.get());
		current = std::move(worklist.front());
		assert(!worklist.front().get());
		worklist.pop_front();

		
    //llvm::errs() << "********* TRACE *********\n";                               ///////////////////////
		//current->annotation.dump();
		//current->graph.to_dot("");
		// if (current->graph.getExtensionFrom() > 0) return;
		

		/*
		current->graph.to_dot("");
		for (auto& cs : current->in_critical_section)
			llvm::errs() << "CS::ipid-" << cs.first << "_d-" << cs.second;
		llvm::errs() << "\n";
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

		if (!current->in_critical_section.empty()) {
      assert(current->in_critical_section.size() == 1);
			int cs_pid = current->in_critical_section.begin()->first;
			for (auto it = nodesToMutate.begin(); it != nodesToMutate.end(); ) {
				const Node * nd = *it;
				if (nd->getEvent()->iid.get_pid() == cs_pid)
					++it;
				else
					it = nodesToMutate.erase(it);
			}
			assert(nodesToMutate.size() == 1);
		}

		// ordering of nodes to try mutations
		// in this branch we DON'T HAVE the preference node
		auto orderedNodesToMutate = std::list<const Node *>();
    for (auto& ndtomut : nodesToMutate) {
      if (ndtomut->getProcessID() == current->graph.starRoot())
				orderedNodesToMutate.push_front(ndtomut); // root first
			else
				orderedNodesToMutate.push_back(ndtomut); // non-root after
		}
		
	  std::vector<unsigned> processLengths = current->graph.getProcessLengths();

		assert(current->unannot.empty());
		for (auto& nd : nodesToMutate)
			current->unannot.insert(nd->getEvent()->iid.get_pid());

    clock_t init = std::clock();
		
		// Get partial-order refinements that order extension events
		// Each refinement will be a candidate for possible mutations
    std::list<PartialOrder> extendedPOs = extensionWritesOrderings();
		assert(!extendedPOs.empty() && "current->trace is one witness");

    time_maz += (double)(clock() - init)/CLOCKS_PER_SEC;
		
    while (!extendedPOs.empty()) {
      auto po = std::move(extendedPOs.front());
			assert(!extendedPOs.front().first.get() &&
						 !extendedPOs.front().second.get());
			extendedPOs.pop_front();


			init = std::clock();
			
			auto withoutMutation = VCValClosure(current->graph, current->annotation);
			//llvm::errs() << "\nclosure after extension writes orderings...\n";
			withoutMutation.valClose(po, nullptr, nullptr);

			time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;
      if (!withoutMutation.closed) ++cl_ordering_failed;
			else ++cl_ordering_succeeded;
			
			if (!withoutMutation.closed) {
        // This ordering of extension events
				// makes the original annotation unrealizable
				// therefore no need to try any mutations
				//llvm::errs() << "Invalidates annotation we have so far\n";                 // ///////////////////////
				po.first.reset();
				po.second.reset();
				continue;
			}
			
			//llvm::errs() << "********* EXTENSION *********\n";                           ///////////////////////
			//current->graph.to_dot(po, "");

			// Try all possible mutations
			for (auto ndit = orderedNodesToMutate.begin();
					 ndit != orderedNodesToMutate.end(); ++ndit) {
				const Node * nd = *ndit;
        if (isRead(nd->getEvent())) {
          bool error = mutateRead(po, withoutMutation, nd);
					if (error) {
            assert(originalTB.error_trace);
						current.reset();
						worklist.clear();
						return error;
					}
				  current->negative.update(nd, processLengths);
				}
				else {
          assert(isLock(nd->getEvent()));
					bool error = mutateLock(po, withoutMutation, nd);
					if (error) {
            assert(originalTB.error_trace);
						current.reset();
						worklist.clear();
						return error;
					}
				}
			}

			po.first.reset();
			po.second.reset();
		}			

		// Delete managed VCTrace
		current.reset();
	}

	return false;
}

/* *************************** */
/* EXTENSION EVENTS ORDERINGS  */
/* *************************** */

std::list<PartialOrder> VCExplorer::extensionWritesOrderings()
{
  assert(current.get());

	current->graph.initWorklist();
	
	// Order extension reads + writes of non-star processes
	// with conflicting writes of non-star processes
	for (unsigned trace_idx = current->graph.getExtensionFrom();
			 trace_idx < current->trace.size(); ++trace_idx) {
		
    const VCEvent& ev = current->trace[trace_idx];
		const Node *nd = current->graph.getNode(ev);
		assert(!isRead(ev) || !current->annotation.defines(nd));
		if (isWrite(ev) && nd->getProcessID() != current->graph.starRoot()) {
			current->graph.orderEventMaz(&ev, current->annotation);
		}
		
	}

	return current->graph.dumpDoneWorklist();
}

std::list<PartialOrder> VCExplorer::readToBeMutatedOrderings(const PartialOrder& po, const Node * nd)
{
  assert(isRead(nd->getEvent()));
  assert(current.get());
	assert(!current->annotation.defines(nd));

	current->graph.initWorklist(po);

	if (nd->getProcessID() != current->graph.starRoot())
	  current->graph.orderEventMaz(nd->getEvent(), current->annotation);

	return current->graph.dumpDoneWorklist();
}

/* *************************** */
/* MUTATE READ                 */
/* *************************** */

bool VCExplorer::mutateRead(const PartialOrder& po, const VCValClosure& withoutMutation, const Node *nd)
{
  assert(isRead(nd->getEvent()));
	
	std::unordered_set<int> mutatedUnannot(current->unannot);
	assert(mutatedUnannot.count(nd->getEvent()->iid.get_pid()));
	mutatedUnannot.erase(nd->getEvent()->iid.get_pid());

  clock_t init = std::clock();
	
  std::list<PartialOrder> readOrderedPOs = readToBeMutatedOrderings(po, nd);
	assert(!readOrderedPOs.empty());

  time_maz += (double)(clock() - init)/CLOCKS_PER_SEC;
	
	while (!readOrderedPOs.empty()) {
		auto roPo = std::move(readOrderedPOs.front());
		assert(!readOrderedPOs.front().first.get() &&
					 !readOrderedPOs.front().second.get());
		readOrderedPOs.pop_front();

    ///////////////////////////////////////////////////////////////////////////////////////////////// SHOULD DO CLOSURE HERE
		
		auto mutationCandidates =
			current->graph.getMutationCandidates(roPo, current->negative, nd);
		
		/*
		llvm::errs() << "NODE[" << nd->getProcessID() << "][" << nd->getEventID() << "] should see:";
		for (auto& valpos_ann : mutationCandidates) {
			llvm::errs() << " " << valpos_ann.first.first << "-";
			char loc = (valpos_ann.first.second == VCAnnotation::Loc::LOCAL)?'L':((valpos_ann.first.second == VCAnnotation::Loc::REMOTE)?'R':'A');
			llvm::errs() << loc;
		} llvm::errs() << "\n\n";
		*/

		for (auto& valpos_ann : mutationCandidates) {

			/*
			llvm::errs() << valpos_ann.first.first << "-";
			char loc = (valpos_ann.first.second == VCAnnotation::Loc::LOCAL)?'L':((valpos_ann.first.second == VCAnnotation::Loc::REMOTE)?'R':'A');
			llvm::errs() << loc << "... ";
			*/

			auto mutatedPo = PartialOrder(std::unique_ptr<ThreadPairsVclocks>
																		(new ThreadPairsVclocks(*(roPo.first))),
																		std::unique_ptr<ThreadPairsVclocks>
																		(new ThreadPairsVclocks(*(roPo.second))));

      init = std::clock();
			
			auto withMutation = VCValClosure(withoutMutation);
			//llvm::errs() << "\nclosure after read orderings and mutation...\n";
			withMutation.valClose(mutatedPo, nd, &(valpos_ann.second));

			time_closure += (double)(clock() - init)/CLOCKS_PER_SEC;
      if (!withMutation.closed) ++cl_mutation_failed;
			else ++cl_mutation_succeeded;
			
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
			mutatedAnnotation.add(nd, valpos_ann.second);

      init = std::clock();
			
			auto mutatedTrace =
				extendTrace(current->graph.linearize(mutatedPo, nullptr),
										mutatedUnannot);

			time_replaying += (double)(clock() - init)/CLOCKS_PER_SEC;
			
			if (mutatedTrace.first.empty())
				return true; // found an error
			assert(traceRespectsAnnotation(mutatedTrace.first, mutatedAnnotation));

      init = std::clock();
			
			VCGraphVclock mutatedGraph(current->graph,       // base for graph
																 std::move(mutatedPo), // base for 'original' po
																 mutatedTrace.first); // to extend the graph
			assert(!mutatedPo.first.get() && !mutatedPo.second.get());

			std::unique_ptr<VCTrace> mutatedVCTrace
				(new VCTrace(std::move(mutatedTrace.first),
										 std::move(mutatedAnnotation),
										 current->negative,
										 std::move(mutatedGraph),
										 std::move(mutatedTrace.second)));	
			assert(mutatedTrace.first.empty() && mutatedAnnotation.empty()
						 && mutatedGraph.empty() && mutatedTrace.second.empty());

      time_graphcopy += (double)(clock() - init)/CLOCKS_PER_SEC;
			
			worklist.push_front(std::move(mutatedVCTrace));
			assert(!mutatedVCTrace.get());
		}
	}

	return false;	
}

/* *************************** */
/* MUTATE LOCK                 */
/* *************************** */

bool VCExplorer::mutateLock(const PartialOrder& po, const VCValClosure& withoutMutation, const Node *nd)
{
  assert(isLock(nd->getEvent()));

  //current->graph.to_dot("");
	
  auto lastLock = current->annotation.getLastLock(nd);

	//llvm::errs() << "\n\nMutating_lock_" << nd->getEvent()->ml.addr.to_string() << ".... ";
	
  if (!lastLock.first) {
    // This lock hasn't been touched before
		// Trivially realizable

		//llvm::errs() << "TRIVIAL!\n\n";
		
		std::unordered_set<int> mutatedUnannot(current->unannot);
		assert(mutatedUnannot.count(nd->getEvent()->iid.get_pid()));
		mutatedUnannot.erase(nd->getEvent()->iid.get_pid());
		
		auto mutatedPo = PartialOrder(std::unique_ptr<ThreadPairsVclocks>
																	(new ThreadPairsVclocks(*(po.first))),
																	std::unique_ptr<ThreadPairsVclocks>
																	(new ThreadPairsVclocks(*(po.second))));
		
		VCAnnotation mutatedAnnotation(current->annotation);
		mutatedAnnotation.setLastLock(nd);

    // current->graph.to_dot(mutatedPo, "");

		auto mutatedLock = VCIID(nd->getProcessID(), nd->getEventID());
		auto mutatedTrace =
			extendTrace(current->graph.linearize(mutatedPo, &mutatedLock),
									mutatedUnannot);
		if (mutatedTrace.first.empty())
			return true; // found an error
		assert(traceRespectsAnnotation(mutatedTrace.first, mutatedAnnotation));

		VCGraphVclock mutatedGraph(current->graph,       // base for graph
															 std::move(mutatedPo), // base for 'original' po
														   mutatedTrace.first); // to extend the graph
    assert(!mutatedPo.first.get() && !mutatedPo.second.get());

		std::unique_ptr<VCTrace> mutatedVCTrace
			(new VCTrace(std::move(mutatedTrace.first),
									 std::move(mutatedAnnotation),
									 current->negative,
									 std::move(mutatedGraph),
									 std::move(mutatedTrace.second)));
		assert(mutatedTrace.first.empty() && mutatedAnnotation.empty()
					 && mutatedGraph.empty() && mutatedTrace.second.empty());

		worklist.push_front(std::move(mutatedVCTrace));
		assert(!mutatedVCTrace.get());
		return false;
	}

	

	assert(lastLock.first);
	auto lastunlockit = current->graph.nodes_iterator(lastLock.second);
	#ifndef NDEBUG
  const Node * lastlocknd = *lastunlockit;
	//llvm::errs() << "lastlock->[ " << lastlocknd->getProcessID() << " ][ " << lastlocknd->getEventID() << " ] ..... ";
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
		//llvm::errs() << "trivfailed!\n";
		return false;
	}
	
	assert(lastunlocknd);

	//llvm::errs() << "lastunlock->[ " << lastunlocknd->getProcessID() << " ][ " << lastunlocknd->getEventID() << " ] ..... ";

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
		//llvm::dbgs() << "failed!\n";
		return false;
	}

	// The mutation on 'mutatedPo' succeeded

	//llvm::errs() << "success!! LASTLOCK OF " << nd->getEvent()->ml.addr.to_string() << " WILL NOW BE [ "
	//						 << nd->getProcessID() << " ][ " << nd->getEventID() << " ]\n\n";

	std::unordered_set<int> mutatedUnannot(current->unannot);
	assert(mutatedUnannot.count(nd->getEvent()->iid.get_pid()));
	mutatedUnannot.erase(nd->getEvent()->iid.get_pid());
	
	VCAnnotation mutatedAnnotation(current->annotation);
	mutatedAnnotation.setLastLock(nd);

  //current->graph.to_dot(mutatedPo, "");
	
	auto mutatedLock = VCIID(nd->getProcessID(), nd->getEventID());
	auto mutatedTrace =
		extendTrace(current->graph.linearize(mutatedPo, &mutatedLock),
								mutatedUnannot);
	if (mutatedTrace.first.empty())
		return true; // found an error
	assert(traceRespectsAnnotation(mutatedTrace.first, mutatedAnnotation));

	VCGraphVclock mutatedGraph(current->graph,       // base for graph
														 std::move(mutatedPo), // base for 'original' po
														 mutatedTrace.first); // to extend the graph
	assert(!mutatedPo.first.get() && !mutatedPo.second.get());

	std::unique_ptr<VCTrace> mutatedVCTrace
		(new VCTrace(std::move(mutatedTrace.first),
								 std::move(mutatedAnnotation),
								 current->negative,
								 std::move(mutatedGraph),
								 std::move(mutatedTrace.second)));	
	assert(mutatedTrace.first.empty() && mutatedAnnotation.empty()
				 && mutatedGraph.empty() && mutatedTrace.second.empty());

	worklist.push_front(std::move(mutatedVCTrace));
	assert(!mutatedVCTrace.get());

	return false;
}

/* *************************** */
/* EXTEND TRACE                */
/* *************************** */

std::pair<std::vector<VCEvent>,
					std::unordered_map<int, int>>
VCExplorer::extendTrace(std::vector<VCEvent>&& tr,
												const std::unordered_set<int>& unannot)
{
	VCTraceBuilder TB(originalTB.config, originalTB.M, std::move(tr), unannot);
	auto trace_cs = TB.extendGivenTrace();

	if (TB.has_error()) {
		// ERROR FOUND
    originalTB.error_trace = TB.get_trace();
		return {std::vector<VCEvent>(), std::unordered_map<int, int>()};
	}
	
	++executed_traces;

	return trace_cs;
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
