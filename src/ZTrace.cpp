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


/* *************************** */
/* TRACE EXTENSION             */
/* *************************** */

ZTraceExtension::ZTraceExtension
(std::vector<ZEvent>&& extension, int ext_from_id,
 bool has_assume_blocked_thread, bool has_deadlocked_thread)
 : extension(std::move(extension)),
   ext_from_id(ext_from_id),
   has_assume_blocked_thread(has_assume_blocked_thread),
   has_deadlocked_thread(has_deadlocked_thread)
{
  assert(extension.empty());
}


std::string ZTraceExtension::to_string() const
{
  std::stringstream res;

  res << trace_to_string(extension);
  res << "Extension starts from id ::: " << ext_from_id << "\n";
  res << "Assume-blocked thread ::: " << has_assume_blocked_thread << "\n"
      << "   Deadlocked thread ::: " << has_deadlocked_thread << "\n";

  return res.str();
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZTraceExtension& ext)
{
  out << ext.to_string();
  return out;
}


void ZTraceExtension::dump() const
{
  llvm::errs() << *this << "\n";
}


/* *************************** */
/* TRACE                       */
/* *************************** */

ZTrace::ZTrace
(const ZTrace *parent_trace,
 std::vector<ZEvent>&& new_exec,
 std::vector<std::unique_ptr<ZEvent>>&& new_tau,
 ZAnnotation&& new_annotation,
 std::set<ZEventID>&& new_committed,
 const ZEventID& now_mutated)
 : parent(parent_trace),
   exec(std::move(new_exec)),
   tau(std::move(new_tau)),
   annotation(std::move(new_annotation)),
   committed(std::move(new_committed)),
   ext_from_id(-1),
   ext_reads_locks(),
   proc_seq_to_thread_id(),
   previously_mutated(now_mutated)
{
  if (parent_trace)
    proc_seq_to_thread_id = std::unordered_map<
      std::vector<int>, unsigned>(parent_trace->proc_seq_to_thread_id);
  assert(new_exec.empty() && new_tau.empty() &&
         new_annotation.empty() && new_committed.empty());
}


std::string ZTrace::to_string(bool noexec = true) const
{
  std::stringstream res;

/*
  if (!parent) {
    res << "########################\n"
        << "#     INITIAL TRACE    #\n"
        << "########################\n";
  } else if (depth > 0) {
    res << parent->to_string(depth-1);
  }
*/
  if (!noexec) {
    res << "EXECUTION\n" << trace_to_string(exec);
  }
  res << "TAU\n" << trace_to_string(tau);
  res << annotation.to_string();
  res << "COMMITTED\n";
  for (const auto& id : committed)
    res << id.to_string() << "  ";
  res << "Extension from: " << ext_from_id
      << " ::: reads/locks in extension:";
  for (int id : ext_reads_locks)
    res << " " << id;
  res << "\n";
  return res.str();
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZTrace& tr)
{
  out << tr.to_string();
  return out;
}


void ZTrace::dump() const
{
  llvm::errs() << *this << "\n";
}
