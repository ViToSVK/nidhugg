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

#ifndef __VC_TRACE_H__
#define __VC_TRACE_H__

#include <vector>

#include "VCAnnotation.h"
#include "VCGraphVclock.h"

class VCTrace {
 public:

  /* *************************** */
  /* INFORMATION                 */
  /* *************************** */
  
  std::vector<VCEvent> trace;

	VCAnnotation annotation;
	//VCAnnotationNeg annotationNeg;

	VCGraphVclock graph;

	std::unordered_set<int> unannot;
	

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  VCTrace() = delete;
	
  VCTrace(std::vector<VCEvent>&& trace,
				  int star_root_index)
	: trace(std::move(trace)),
		annotation(),
		graph(this->trace, star_root_index),
		unannot() {};

  VCTrace(std::vector<VCEvent>&& trace,
					VCAnnotation&& annotation,
					VCGraphVclock&& graph)
	: trace(std::move(trace)),
		annotation(std::move(annotation)),
		graph(std::move(graph)),
		unannot() {};
	
  VCTrace(VCTrace&& tr) = default;
  VCTrace& operator=(VCTrace&& tr) = delete;
  VCTrace(const VCTrace&) = delete;
  VCTrace& operator=(const VCTrace&) = delete;
  
};

#endif
