/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2020 Viktor Toman
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

#include "ZGraph.h"


class ZTrace {
 private:
  const ZTrace *_parent;

  const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>> _trace;
  const ZAnnotation _annotation;
  ZAnnotationNeg _negative;
  const std::shared_ptr<ZGraph> _graph;
  const std::shared_ptr<ZPartialOrder> _po_part;

 public:
  const std::vector<std::unique_ptr<ZEvent>>& trace() const;
  const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& trace_ptr() const;
  const ZAnnotation& annotation() const;
  ZAnnotationNeg& negative();
  const ZAnnotationNeg& negative() const;
  const ZGraph& graph() const;
  const ZPartialOrder& po_part() const;

  void set_negative(const ZAnnotationNeg& oth);

  // Whether trace has an assume-blocked thread
  bool assumeblocked;
  // Whether this annotated trace is not full but no mutation
  // is possible (i.e. deadlocked). We'll count it as full
  mutable bool deadlocked;

  std::map<const ZEvent *, std::map<int, std::unique_ptr<ZTrace>>> children_read;
  std::map<const ZEvent *, std::unique_ptr<ZTrace>> children_lock;

  bool empty() const;
  std::string to_string(unsigned depth) const;
  void dump() const;

  /* *************************** */
  /* HELPERS                     */
  /* *************************** */

  std::list<const ZEvent *> events_to_mutate() const {
    return graph().events_to_mutate(annotation());
  }

  std::set<ZAnn> mutation_candidates(const ZEvent *read) const {
    return graph().mutation_candidates_grouped(po_part(), read, negative());
  }

  /* *************************** */
  /* CONSTRUCTORS                */
  /* *************************** */

  ZTrace() = delete;

  ZTrace(const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& initial_trace,
         bool assumeblocked);

  ZTrace(const ZTrace& parentTrace,
         const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& new_trace,
         ZAnnotation&& new_annotation,
         ZPartialOrder&& new_po,
         bool assumeblocked);

  ZTrace(ZTrace&& tr) = default;
  ZTrace(const ZTrace& tr) = delete;
  ZTrace& operator=(ZTrace&& tr) = delete;
  ZTrace& operator=(const ZTrace& tr) = delete;

};

#endif // __Z_TRACE_H__
