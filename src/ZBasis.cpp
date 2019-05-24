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


const LineT& ZBasis::operator[](std::pair<unsigned, int> ids) const
{
  auto it = thread_aux_to_line_id.find(ids);
  assert(it != thread_aux_to_line_id.end());
  return lines[it->second];
}


const LineT& ZBasis::operator[](unsigned thread_id, int aux_id) const
{
  return operator[](std::pair<unsigned, int>(thread_id, aux_id));
}


const ZEvent * ZBasis::getEvent(unsigned thread_id, int aux_id, unsigned event_id) const
{
  assert(hasThreadAux(thread_id, aux_id));
  const LinesT& line = lines[thread_id, aux_id];
  assert(event_id < line.size());
  assert(hasEvent(line[event_id]));
  return line[event_id];
}


void ZBasis::addLine(const ZEvent * ev)
{
  assert(!hasEvent(ev));
  assert(ev->thread_id < 20 && "Thread ID not set up yet");
  assert(!hasThreadAux(ev->threadID(), ev->auxID()));

  thread_aux_to_line_id.emplace(std::pair<ev->threadID(), ev->auxID()>, lines.size());
  lines.push_back(LineT());
  lines.back().reserve(8);
}


void ZBasis::addEvent(const ZEvent * ev)
{
  assert(!hasEvent(ev));
  assert(ev->thread_id < 20 && "Thread ID not set up yet");
  assert(hasThreadAux(ev->threadID(), ev->auxID()));

  auto it = thread_aux_to_line_id.find(ev->threadID(), ev->auxID());
  assert(it != thread_aux_to_line_id.end());
  LineT& line = lines[it->second];

  event_id = line.size();
  event_to_position.emplace(ev, std::pair<unsigned,unsigned>(it->second, event_id));
  line.push_back(ev);
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
}


bool ZBasis::hasThreadAux(std::pair<unsigned, int> ids) const
{
  return thread_aux_to_line_id.count(ids);
}


bool ZBasis::hasThreadAux(unsigned thread_id, int aux_id) const
{
  return hasThreadAux(std::pair<unsigned, int>(thread_id, aux_id));
}


std::unordered_map<std::pair<unsigned, int>, unsigned> ZBasis::line_sizes() const
{
  auto res = std::unordered_map<std::pair<unsigned, int>, unsigned>();
  for (const auto& thaux_line : thread_aux_to_line_id) {
    auto thaux = thaux_line.first;
    unsigned size = (*this)[thaux].size();
    res.emplace(thaux, size);
  }
  return res;
}


// <threadID, newly_added?>
std::pair<unsigned, bool> ZBasis::getThreadID(const ZEvent * ev)
{
  auto proc_seq = ev->cpid.get_proc_seq();
  auto it = proc_seq_to_thread_id.find(proc_seq);
  if (it == proc_seq_to_thread_id.end()) {
    unsigned res = proc_seq_to_thread_id.size();
    proc_seq_to_thread_id.emplace_hint(it, proc_seq, res);
    return {res, true};
  };
  return {it->second, false};
}


bool ZBasis::hasEvent(const ZEvent *ev) const
{
  return event_to_position.find(ev) != event_to_position.end();
}
