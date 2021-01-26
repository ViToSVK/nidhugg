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
 const std::shared_ptr<ZGraph>& initial_graph,
 const std::shared_ptr<ZPartialOrder>& initial_po_full,
 const std::shared_ptr<ZPartialOrder>& initial_po_part,
 bool assumeblocked)
  : _parent(nullptr),
    _trace(initial_trace),
    _annotation(),
    _negative(),
    _graph(initial_graph),
    _po_full(initial_po_full),
    _po_part(initial_po_part),
    assumeblocked(assumeblocked),
    deadlocked(false)
{
  assert(_trace && _graph && _po_full && _po_part);
  assert(!empty());
  assert(&(_po_full->graph) == _graph.get());
  assert(&(_po_part->graph) == _graph.get());
}


ZTrace::ZTrace
(const ZTrace& parentTrace,
 const std::shared_ptr<std::vector<std::unique_ptr<ZEvent>>>& new_trace,
 ZAnnotation&& new_annotation,
 const std::shared_ptr<ZGraph>& new_graph,
 const std::shared_ptr<ZPartialOrder>& new_po_full,
 const std::shared_ptr<ZPartialOrder>& new_po_part,
 bool assumeblocked)
  : _parent(&parentTrace),
    _trace(new_trace),
    _annotation(std::move(new_annotation)),
    _negative(parentTrace.negative()),
    _graph(new_graph),
    _po_full(new_po_full),
    _po_part(new_po_part),
    assumeblocked(assumeblocked),
    deadlocked(false)
{
  assert(_trace && _graph && _po_full && _po_part);
  assert(!empty());
  assert(&(_po_full->graph) == _graph.get());
  assert(&(_po_part->graph) == _graph.get());
}


bool ZTrace::empty() const
{
  return ((!_trace || _trace->empty()) &&
          annotation().empty() &&
          _negative.empty() &&
          (!_graph || _graph->empty()) &&
          (!_po_full || _po_full->empty()) &&
          (!_po_part || _po_part->empty()));
}


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


const std::shared_ptr<ZGraph>& ZTrace::graph_ptr() const
{
  assert(_graph);
  return _graph;
}


const ZPartialOrder& ZTrace::po_full() const
{
  assert(_po_full);
  return *_po_full;
}


const std::shared_ptr<ZPartialOrder>& ZTrace::po_full_ptr() const
{
  assert(_po_full);
  return _po_full;
}


const ZPartialOrder& ZTrace::po_part() const
{
  assert(_po_part);
  return *_po_part;
}


void ZTrace::set_negative(const ZAnnotationNeg& oth)
{
  assert(negative().empty());
  _negative.set_mapping(oth);
}


std::string ZTrace::to_string(int depth = 2) const
{
  std::stringstream res;
  assert(depth >= 0);

  if (!_parent) {
    res << "########################\n"
        << "#     INITIAL TRACE    #\n"
        << "########################\n";
  } else if (depth > 0) {
    res << _parent->to_string(depth-1);
  }

  res << po_part().to_string();
  res << annotation().to_string();
  res << "\nvvvvvvvvvvvvvvvvvvvvvvvv\n";

  return res.str();
}


void ZTrace::dump() const
{
  llvm::errs() << to_string() << "\n";
}
