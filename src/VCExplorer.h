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

#ifndef __VC_EXPLORER_H__
#define __VC_EXPLORER_H__

#include <list>
#include <memory>

#include "VCTraceBuilder.h"
#include "VCTrace.h"

class VCExplorer {

  VCTraceBuilder& originalTB;
	
  std::list<std::unique_ptr<VCTrace>> worklist;

	std::unique_ptr<VCTrace> current;

	//

  std::list<PartialOrder> extensionWritesOrderings();
	
  std::vector<VCEvent> extendTrace(std::vector<VCEvent>&& tr,
																	 const std::unordered_set<int>& unannot);

	std::unordered_set<const Node *> getNodesToMutate();

	void mutateRead(const PartialOrder& po, const Node *nd);

	void mutateLock(const PartialOrder& po, const Node *nd);
	
  /* *************************** */
  /* STATISTICS                  */
  /* *************************** */ 
  
  // Number of executed traces (1 because of the initial
  // trace that was obtained as a constructor argument)
  unsigned executed_traces = 1;
  // Total number of executed instructions over all traces
  // TODO: count only over annotation-leaf traces???
  unsigned total_instr_executed = 0;
  
 public:

  void explore();

  void print_stats();
  
  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  VCExplorer(std::vector<VCEvent>&& trace, VCTraceBuilder& tb)
		: originalTB(tb) {
    worklist.push_back(std::unique_ptr<VCTrace>(new VCTrace(std::move(trace))));
  }
  
};

#endif
