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
#include "VCTraceBuilder.h"
#include "VCHappensBeforeGraph.h"


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

		assert(!current.get());
		current = std::move(worklist.front());
		assert(!worklist.front().get());
		worklist.pop_front();

		// Compute the basis
		assert(current->basis.empty());
		current->basis = VCBasis(current->trace);

		// Initialize the HappensBeforeGraph
		VCHappensBeforeGraph hbgraph(current->basis, current->annotation);
		

		// Delete managed VCTrace
		current.reset();
	}
}
