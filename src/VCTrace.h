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

  std::unordered_set<int> unannot;

  std::unordered_map<int, int> in_critical_section;

  unsigned processMutationPreference;

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  VCTrace() = delete;

  VCTrace(std::vector<VCEvent>&& trace,
          int star_root_index)
  : trace(std::move(trace)),
    annotation(),
    negative(),
    graph(this->trace, star_root_index),
    unannot(),
    in_critical_section(),
    processMutationPreference(0)
      {
        for (unsigned i = 0; i < trace.size(); ++i)
          if (isRead(trace[i])) {
            const Node *nd = graph.getNode(trace[i]);
            if (!graph.nodes_iterator(nd).atProcessEnd()) {
              assert(nd->getProcessID() == 0);
              annotation.add(nd, VCAnnotation::Ann());
            }
          }
      };

  VCTrace(std::vector<VCEvent>&& trace,
          VCAnnotation&& annotation,
          const VCAnnotationNeg& negative,
          VCGraphVclock&& graph,
          std::unordered_map<int, int>&& cs,
          unsigned pref)
  : trace(std::move(trace)),
    annotation(std::move(annotation)),
    negative(negative),
    graph(std::move(graph)),
    unannot(),
    in_critical_section(std::move(cs)),
    processMutationPreference(pref)
      {};

  VCTrace(VCTrace&& tr) = default;
  VCTrace& operator=(VCTrace&& tr) = delete;
  VCTrace(const VCTrace&) = delete;
  VCTrace& operator=(const VCTrace&) = delete;

};

#endif
