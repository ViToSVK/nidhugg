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

#include "ZTrace.h"


ZTrace::ZTrace
(const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& initial_trace,
 bool assumeblocked)
  : _parent(nullptr),
    _trace(initial_trace),
    _annotation(),
    _negative(),
    _graph(new ZGraph(this->trace())),
    assumeblocked(assumeblocked),
    deadlocked(false)
{}


ZTrace::ZTrace
(const ZTrace& parentTrace,
 const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& new_trace,
 ZAnnotation&& new_annotation,
 ZPartialOrder&& new_po,
 bool assumeblocked)
  : _parent(&parentTrace),
    _trace(new_trace),
    _annotation(std::move(new_annotation)),
    _negative(parentTrace.negative()),
    _graph(new ZGraph(parentTrace.graph(), std::move(new_po),
                      this->trace(), this->annotation())),
    assumeblocked(assumeblocked),
    deadlocked(false)
{}


const std::vector<std::unique_ptr<ZEvent>>& ZTrace::trace() const
{
  assert(_trace);
  return *_trace;
}


const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& ZTrace::trace_ptr() const
{
  assert(_trace);
  return _trace;
}

const ZAnnotation& ZTrace::annotation() const { return _annotation; }
ZAnnotationNeg& ZTrace::negative() { return _negative; }
const ZAnnotationNeg& ZTrace::negative() const { return _negative; }

const ZGraph& ZTrace::graph() const
{
  assert(_graph);
  return *_graph;
}


bool ZTrace::empty() const
{
  return (!_trace && annotation().empty() &&
          _negative.empty() && !_graph);
}


std::string ZTrace::to_string(unsigned depth = 2) const
{
  std::stringstream res;

  if (!_parent) {
    res << "########################\n"
        << "#     INITIAL TRACE    #\n"
        << "########################\n";
  } else if (depth > 0) {
    res << _parent->to_string(depth-1);
  }

  res << graph().to_string();
  res << annotation().to_string();
  res << "\nvvvvvvvvvvvvvvvvvvvvvvvv\n";

  return res.str();
}


void ZTrace::dump() const
{
  llvm::errs() << to_string() << "\n";
}
