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

#include "ZBasis.h"
#include "ZGraph.h"


ZGraph ZBasis::graphDummy = ZGraph();


// Empty
ZBasis::ZBasis(const ZGraph& gr)
  : graph(gr),
    root_thread_id(INT_MAX),
    init(ZEvent(true)),
    lines(),
    thread_aux_to_line_id(),
    proc_seq_to_thread_id(),
    event_to_position()
{
  assert(graph.empty() && empty());
}


// Empty
ZBasis::ZBasis()
  : ZBasis(ZBasis::graphDummy) {}


// Initial
ZBasis::ZBasis(const ZGraph& gr, int root_thread_id)
  : graph(gr),
    root_thread_id(root_thread_id),
    init(ZEvent(true)),
    lines(),
    thread_aux_to_line_id(),
    proc_seq_to_thread_id(),
    event_to_position()
{
  assert(empty());
}


// When extending
ZBasis::ZBasis(const ZBasis& oth, const ZGraph& gr)
  : graph(gr),
    root_thread_id(oth.root_thread_id),
    init(ZEvent(true)),
    lines(oth.lines),
    thread_aux_to_line_id(oth.thread_aux_to_line_id),
    threads_auxes(oth.threads_auxes),
    proc_seq_to_thread_id(oth.proc_seq_to_thread_id),
    event_to_position(oth.event_to_position)
{}


const LineT& ZBasis::operator()(std::pair<unsigned, int> ids) const
{
  assert(ids.first != INT_MAX && "Called for initial event");
  auto it = thread_aux_to_line_id.find(ids);
  assert(it != thread_aux_to_line_id.end());
  return lines[it->second];
}


const LineT& ZBasis::operator()(unsigned thread_id, int aux_id) const
{
  return operator()(std::pair<unsigned, int>(thread_id, aux_id));
}


const ZEvent * ZBasis::getEvent(unsigned thread_id, int aux_id, unsigned event_id) const
{
  assert(thread_id != INT_MAX && "Called for initial event");
  assert(hasThreadAux(thread_id, aux_id));
  const LineT& line = this->operator()(thread_id, aux_id);
  assert(event_id < line.size());
  assert(hasEvent(line[event_id]));
  return line[event_id];
}


const ZEvent * ZBasis::getEvent(const ZObs& obs) const
{
  return getEvent(obs.thr, -1, obs.ev);
}


const ZEvent * ZBasis::getUnlockOfThisLock(const ZObs& obs) const
{
  unsigned curEv = obs.ev;
  auto lock = getEvent(obs.thr, -1, curEv);
  assert(isLock(lock));

  while (curEv + 1 < (*this)(obs.thr, -1).size()) {
    curEv++;
    auto res = getEvent(obs.thr, -1, curEv);
    if (isUnlock(res) && sameMl(lock, res))
      return res;
  }

  return nullptr;
}


void ZBasis::addLine(const ZEvent *ev)
{
  assert(!hasEvent(ev));
  assert(ev->threadID() < 20 && "Thread ID not set up yet");
  assert(!hasThreadAux(ev->threadID(), ev->auxID()));

  auto key = std::pair<unsigned, int>(ev->threadID(), ev->auxID());
  thread_aux_to_line_id.emplace(key, lines.size());
  assert(ev->threadID() <= threads_auxes.size());
  if (ev->threadID() == threads_auxes.size())
    threads_auxes.push_back(std::set<int>());
  assert(!threads_auxes.at(ev->threadID()).count(ev->auxID()));
  threads_auxes[ev->threadID()].emplace(ev->auxID());
  lines.push_back(LineT());
  lines.back().reserve(16);
}


void ZBasis::addEvent(const ZEvent *ev)
{
  assert(!hasEvent(ev));
  assert(ev->threadID() < 20 && "Thread ID not set up yet");
  assert(hasThreadAux(ev->threadID(), ev->auxID()));

  auto key = std::pair<unsigned, int>(ev->threadID(), ev->auxID());
  auto it = thread_aux_to_line_id.find(key);
  assert(it != thread_aux_to_line_id.end());
  LineT& line = lines[it->second];

  unsigned event_id = line.size();
  assert(event_id == ev->eventID());
  event_to_position.emplace(ev, std::pair<unsigned,unsigned>(it->second, event_id));
  line.push_back(ev);
  assert(hasEvent(ev));
}


void ZBasis::replaceEvent(const ZEvent *oldEv, const ZEvent *newEv)
{
  assert(oldEv != newEv && (*oldEv == *newEv)
         && "different pointers but same thread/aux/eventID");
  assert(oldEv->kind == newEv->kind &&
         oldEv->cpid == newEv->cpid &&
         oldEv->ml == newEv->ml);
  assert(hasEvent(oldEv) && !hasEvent(newEv));

  auto line_event = event_to_position[oldEv];
  event_to_position.erase(oldEv);
  event_to_position.emplace(newEv, line_event);

  assert(line_event.first < lines.size() &&
         line_event.second < lines[line_event.first].size() &&
         lines[line_event.first][line_event.second] == oldEv);
  lines[line_event.first][line_event.second] = newEv;
  assert(!hasEvent(oldEv) && hasEvent(newEv));
}


void ZBasis::shrink()
{
  lines.shrink_to_fit();
  for (unsigned i=0; i<lines.size(); ++i)
    lines[i].shrink_to_fit();
}


unsigned ZBasis::lineID(unsigned thread_id, int aux_id) const
{
  assert(thread_id != INT_MAX && "Called for initial event");
  assert(hasThreadAux(thread_id, aux_id));
  auto key = std::pair<unsigned,int>(thread_id, aux_id);
  auto it = thread_aux_to_line_id.find(key);
  assert(it != thread_aux_to_line_id.end());
  return it->second;
}


unsigned ZBasis::lineID(const ZEvent *ev) const
{
  return lineID(ev->threadID(), ev->auxID());
}


bool ZBasis::hasThreadAux(std::pair<unsigned, int> ids) const
{
  assert(ids.first != INT_MAX && "Called for initial event");
  assert(threads_auxes.size() >= ids.first ||
         !threads_auxes.at(ids.first).count(ids.second) ||
         thread_aux_to_line_id.count(ids));
  assert(!thread_aux_to_line_id.count(ids) ||
         threads_auxes.at(ids.first).count(ids.second));
  return thread_aux_to_line_id.count(ids);
}


bool ZBasis::hasThreadAux(unsigned thread_id, int aux_id) const
{
  return hasThreadAux(std::pair<unsigned, int>(thread_id, aux_id));
}


std::vector<unsigned> ZBasis::real_sizes_minus_one() const
{
  std::vector<unsigned> res;
  for (unsigned thr = 0; thr < number_of_threads(); ++thr) {
    assert(!(*this)(thr, -1).empty());
    res.push_back((*this)(thr, -1).size() - 1);
  }
  return res;
}


unsigned ZBasis::number_of_threads() const
{
  return threads_auxes.size();
}


const std::set<int>& ZBasis::auxes(unsigned thread_id) const
{
  assert(thread_id != INT_MAX && "Called for initial event");
  assert(thread_id < threads_auxes.size());
  return threads_auxes[thread_id];
}


int ZBasis::auxForMl(const SymAddrSize& ml, unsigned thr) const
{
  assert(thr < number_of_threads());
  auto axs = auxes(thr);
  assert(!axs.empty());
  if (axs.size() == 1) {
    // No write events in this thread
    return -1;
  }
  if (axs.size() == 2 && graph.tso) {
    #ifndef NDEBUG
    auto it = axs.begin();
    while (*it == -1) {
      ++it;
      assert(it != axs.end());
    }
    assert(*it == 0);
    #endif
    return 0;
  }
  // PSO below
  for (auto aux : axs) {
    assert(!((*this)(thr, aux).empty()));
    const ZEvent *first = (*this)(thr, aux)[0];
    if (isWriteM(first) && first->ml == ml)
      return aux;
  }
  // This thread has no event for this ml
  return -1;
}


// <threadID, added_with_this_call?>
std::pair<unsigned, bool> ZBasis::getThreadID(const std::vector<int>& proc_seq)
{
  auto it = proc_seq_to_thread_id.find(proc_seq);
  if (it == proc_seq_to_thread_id.end()) {
    unsigned res = proc_seq_to_thread_id.size();
    proc_seq_to_thread_id.emplace_hint(it, proc_seq, res);
    return {res, true};
  };
  return {it->second, false};
}


// <threadID, added_with_this_call?>
std::pair<unsigned, bool> ZBasis::getThreadID(const ZEvent *ev)
{
  assert(!isInitial(ev) && "Called for initial event");
  return getThreadID(ev->cpid.get_proc_seq());
}


unsigned ZBasis::getThreadIDnoAdd(const std::vector<int>& proc_seq) const
{
  auto it = proc_seq_to_thread_id.find(proc_seq);
  if (it == proc_seq_to_thread_id.end())
    return 1337;
  else
    return it->second;
}


unsigned ZBasis::getThreadIDnoAdd(const ZEvent * ev) const
{
  assert(!isInitial(ev) && "Called for initial event");
  return getThreadIDnoAdd(ev->cpid.get_proc_seq());
}


bool ZBasis::hasEvent(const ZEvent *ev) const
{
  if (ev == initial())
    return true;
  auto it = event_to_position.find(ev);
  if (it == event_to_position.end())
    return false;
  assert(lines[it->second.first][it->second.second] == ev);
  assert(lineID(ev) == it->second.first);
  assert(ev->eventID() == it->second.second);
  return true;
}
