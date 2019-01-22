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

#include "VCAnnotationNeg.h"
#include "VCGraphVclock.h"

class VCTrace {
 public:

  std::vector<VCEvent> trace;

  VCAnnotation annotation;
  VCAnnotationNeg negative;

  VCGraphVclock graph;

  unsigned processMutationPreference;

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  VCTrace() = delete;

  VCTrace(std::vector<VCEvent>&& initial_trace,
          int star_root_index)
  : trace(std::move(initial_trace)),
    annotation(),
    negative(),
    graph(this->trace, star_root_index),
    processMutationPreference(0)
      {};

  VCTrace(std::vector<VCEvent>&& trace,
          const VCAnnotation& annotation,
          const VCAnnotationNeg& negative,
          VCGraphVclock&& graph,
          unsigned pref)
  : trace(std::move(trace)),
    annotation(annotation),
    negative(negative),
    graph(std::move(graph)),
    processMutationPreference(pref)
      {};

  VCTrace(VCTrace&& tr) = default;
  VCTrace& operator=(VCTrace&& tr) = delete;
  VCTrace(const VCTrace&) = delete;
  VCTrace& operator=(const VCTrace&) = delete;

};

#endif
