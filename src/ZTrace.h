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

#include <vector>

#include "ZAnnotation.h"


class ZTraceExtension {
 public:
  std::vector<ZEvent> extension;
  int ext_from_id;
  bool has_assume_blocked_thread;
  bool has_deadlocked_thread;

  ZTraceExtension() = default;
  ZTraceExtension(std::vector<ZEvent>&& extension,
                  int ext_from_id,
                  bool has_assume_blocked_thread,
                  bool has_deadlocked_thread);

  ZTraceExtension(const ZTraceExtension&) = delete;
  ZTraceExtension(ZTraceExtension&&) = default;
  ZTraceExtension& operator=(const ZTraceExtension&) = delete;
  ZTraceExtension& operator=(ZTraceExtension&&) = default;

  bool empty() const { return (extension.empty()); }

  std::string to_string() const;
  void dump() const;
};
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZTraceExtension& ext);


class ZTrace {
 private:
  const ZTrace *parent;

 public:
  std::vector<ZEvent> exec;
  std::vector<std::unique_ptr<ZEvent>> tau;
  ZAnnotation annotation;
  const std::set<ZEventID> committed;
  int ext_from_id;
  std::vector<int> ext_reads_locks;

  ZTrace() = delete;
  ZTrace(const ZTrace *parent_trace,
         std::vector<ZEvent>&& new_exec,
         std::vector<std::unique_ptr<ZEvent>>&& new_tau,
         ZAnnotation&& new_annotation,
         std::set<ZEventID>&& new_commited);

  ZTrace(const ZTrace&) = delete;
  ZTrace(ZTrace&&) = default;
  ZTrace& operator=(const ZTrace&) = delete;
  ZTrace& operator=(ZTrace&&) = delete;

  bool empty() const {
    return (exec.empty() && tau.empty() &&
            annotation.empty() && committed.empty());
  }

  std::string to_string(unsigned depth) const;
  void dump() const;
};
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZTrace& tr);

#endif // __Z_TRACE_H__
