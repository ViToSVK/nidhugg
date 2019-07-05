/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2019 Viktor Toman
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
 private:
  const ZTrace *parent;

 public:

  const std::vector<ZEvent> trace;

  const ZAnnotation annotation;
  ZAnnotationNeg negative;

  ZGraph graph;

  // Whether trace has an assume-blocked thread
  bool assumeblocked;
  // Whether this annotated trace is not full but no mutation
  // is possible (i.e. deadlocked). We'll count it as full
  mutable bool deadlocked;

  bool empty() const {
    return (trace.empty() && annotation.empty() &&
            negative.empty() && graph.empty());
  }

  std::string to_string(unsigned depth) const;
  void dump() const;


  /* *************************** */
  /* HELPERS                     */
  /* *************************** */

  std::list<const ZEvent *> getEventsToMutate() const {
    return graph.getEventsToMutate(annotation);
  }

  std::list<ZObs> getObsCandidates(const ZEvent *read) const {
    return graph.getObsCandidates(read, negative);
  }

  const ZEvent *getEvent(const ZObs& obs) const {
    return graph.getBasis().getEvent(obs);
  }

  bool isRoot(unsigned thr_id) const {
    return graph.getBasis().isRoot(thr_id);
  }

  bool isRoot(const ZEvent *ev) const {
    return graph.getBasis().isRoot(ev);
  }


  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  ZTrace();

  ZTrace(std::vector<ZEvent>&& initial_trace,
         int star_root_index, bool assumeblocked, bool tso);

  ZTrace(const ZTrace& parentTrace,
         std::vector<ZEvent>&& new_trace,
         ZAnnotation&& new_annotation,
         ZPartialOrder&& new_po,
         bool assumeblocked);

  ZTrace(ZTrace&& tr) = default;
  ZTrace& operator=(ZTrace&& tr) = delete;
  ZTrace(const ZTrace&) = delete;
  ZTrace& operator=(const ZTrace&) = delete;

};

#endif // __Z_TRACE_H__
