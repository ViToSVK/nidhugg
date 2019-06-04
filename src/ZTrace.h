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

#ifndef __Z_TRACE_H__
#define __Z_TRACE_H__

#include <vector>

#include "ZAnnotationNeg.h"
#include "ZGraph.h"


class ZTrace {
 public:

  const std::vector<ZEvent> trace;

  const ZAnnotation annotation;
  const ZAnnotationNeg negative;

  ZGraph graph;

  // Whether this trace is not full but no mutation is
  // possible (i.e. deadlocked). We'll count it as full
  mutable bool deadlocked;

  bool empty() const {
    return (trace.empty() && annotation.empty() &&
            negative.empty() && graph.empty());
  }

  std::list<const ZEvent *> getEventsToMutate() const {
    return graph.getEventsToMutate(annotation);
  }

  std::list<ZObs> getObsCandidates(const ZEvent *read) const {
    return graph.getObsCandidates(read, negative);
  }

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  ZTrace()
    : trace(), annotation(), negative(), graph(), deadlocked(false)
      {
        assert(empty());
      };

  ZTrace(std::vector<ZEvent>&& initial_trace,
         int star_root_index)
  : trace(std::move(initial_trace)),
    annotation(),
    negative(),
    graph(this->trace, star_root_index),
    deadlocked(false)
      {};

  ZTrace(std::vector<ZEvent>&& new_trace,
         const ZAnnotation& new_annotation,
         const ZAnnotationNeg& new_negative,
         const ZGraph& old_graph,
         ZPartialOrder&& new_po)
  : trace(std::move(new_trace)),
    annotation(new_annotation),
    negative(new_negative),
    graph(old_graph, std::move(new_po),
          this->trace, this->annotation),
    deadlocked(false)
      {};

  ZTrace(ZTrace&& tr) = default;
  ZTrace& operator=(ZTrace&& tr) = delete;
  ZTrace(const ZTrace&) = delete;
  ZTrace& operator=(const ZTrace&) = delete;

};

#endif // __Z_TRACE_H__
