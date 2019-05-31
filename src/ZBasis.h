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

#ifndef _Z_BASIS_H_
#define _Z_BASIS_H_

#include "ZEvent.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>


namespace std {
  template <typename A, typename B>
  struct hash<std::pair<A, B>> {
    std::size_t operator()(const std::pair<A, B>& pair) const {
      return (hash<A>()(pair.first) << 12) +
             hash<B>()(pair.second);
    }
  };
}


class ZGraph;

typedef std::vector<const ZEvent *> LineT;
typedef std::vector<LineT> LinesT;


class ZBasis {
  friend class ZPartialOrder;
 public:
  // ZGRAPH REFERENCE (given at constructor-time)
  const ZGraph& graph;  ////


  // ROOT (retained accross recursion children)
 private:
  const unsigned root_thread_id; ////
 public:
  unsigned root() const { return root_thread_id; }
  bool isRoot(unsigned idx) const { return root_thread_id == idx; }


  // LINES (changed at recursion child with new ZEvent pointers from new trace)
 private:
  const ZEvent init;
  LinesT lines; ////
 public:
  const ZEvent *initial() const { return &init; }
  const LineT& operator[](std::pair<unsigned, int> ids) const;
  const LineT& operator()(unsigned thread_id, int aux_id) const;
  //
  const ZEvent *getEvent(unsigned thread_id, int aux_id, unsigned event_id) const;
  void addLine(const ZEvent *ev);
  void addEvent(const ZEvent *ev);
  void replaceEvent(const ZEvent *oldEv, const ZEvent *newEv);


  // THREAD_AUX->LINE_ID (retained accross recursion children)
 private:
  std::unordered_map<std::pair<unsigned, int>, unsigned> thread_aux_to_line_id; ////
  unsigned lineID(unsigned thread_id, int aux_id) const;
  unsigned lineID(const ZEvent *ev) const;
 public:
  bool hasThreadAux(std::pair<unsigned, int> ids) const;
  bool hasThreadAux(unsigned thread_id, int aux_id) const;

  std::unordered_map<std::pair<unsigned, int>, unsigned> line_sizes() const;


  // PROCSEQ->THREAD_ID (retained accross recursion children)
 private:
  std::unordered_map<std::vector<int>, unsigned> proc_seq_to_thread_id; ////
 public:
  // <threadID, added_with_this_call?>
  std::pair<unsigned, bool> getThreadID(const ZEvent * ev);


  // EVENT->POSITION (changed at recursion child with new ZEvent pointers from new trace)
 private:
  std::unordered_map<const ZEvent *, std::pair<unsigned, unsigned>> event_to_position; ////
 public:
  bool hasEvent(const ZEvent *ev) const;


 public:
  // Empty
  ZBasis(const ZGraph& gr);
  // Initial
  ZBasis(const ZGraph& gr, int root_thread_id);
  // When extending
  ZBasis(const ZBasis& oth, const ZGraph& gr);

  ZBasis(ZBasis&& oth) = default;
  ZBasis& operator=(ZBasis&& oth) = delete;
  ZBasis(const ZBasis& oth) = delete;
  ZBasis& operator=(const ZBasis& oth) = delete;

  bool empty() const { return lines.empty(); }
  size_t size() const { return lines.size(); }
  size_t events_size() const {
    size_t res = 0;
    for (auto& ln : lines)
      res += ln.size();
    return res;
  }

  // typename clarifies that (const_)iterator
  // is a class and not a static member
  typename LinesT::iterator begin() { return lines.begin(); }
  typename LinesT::iterator end() { return lines.end(); }
  typename LinesT::const_iterator begin() const { return lines.begin(); }
  typename LinesT::const_iterator end() const { return lines.end(); }

};


#endif // _Z_BASIS_H_
