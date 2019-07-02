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

#include "ZTrace.h"


ZTrace::ZTrace()
  : parent(nullptr), trace(), annotation(), negative(),
    graph(), assumeblocked(false), deadlocked(false)
{ assert(empty()); }


ZTrace::ZTrace
(std::vector<ZEvent>&& initial_trace,
 int star_root_index, bool assumeblocked, bool tso)
  : parent(nullptr),
    trace(std::move(initial_trace)),
    annotation(),
    negative(),
    graph(this->trace, star_root_index, tso),
    assumeblocked(assumeblocked),
    deadlocked(false)
{}


ZTrace::ZTrace
(const ZTrace& parentTrace,
 std::vector<ZEvent>&& new_trace,
 ZAnnotation&& new_annotation,
 ZPartialOrder&& new_po,
 bool assumeblocked)
  : parent(&parentTrace),
    trace(std::move(new_trace)),
    annotation(std::move(new_annotation)),
    negative(parentTrace.negative),
    graph(parentTrace.graph, std::move(new_po),
          this->trace, this->annotation),
    assumeblocked(assumeblocked),
    deadlocked(false)
{}


std::string ZTrace::to_string(unsigned depth = 2) const
{
  std::stringstream res;

  if (!parent) {
    res << "########################\n"
        << "#     INITIAL TRACE    #\n"
        << "########################\n";
  } else if (depth > 0) {
    res << parent->to_string(depth-1);
  }

  res << graph.getPo().to_string();
  res << annotation.to_string();
  res << "\nvvvvvvvvvvvvvvvvvvvvvvvv\n\n";

  return res.str();
}


void ZTrace::dump() const
{
  llvm::errs() << to_string();
}
