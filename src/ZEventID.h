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

#ifndef _Z_EVENTID_H_
#define _Z_EVENTID_H_

#include <sstream>

#include "CPid.h"


/* Event identifier, consisting of thread identifier and position within the
 * thread. With a fixed observation function resp. value function, this
 * identifier information uniquely identifies an event.
 * Be careful in PSO, as CPid-aux number may change for a fixed variable-buffer
 * across exploration (e.g. in a fixed thread, variable x was the first to
 * enqueue, hence the respective buffer has aux 0; later in the exploration
 * algorithm, due to annotation change, x is second to enqueue, hence the
 * respective buffer has now aux 1).
 */
class ZEventID {

 public:
  ZEventID() = default;
  // Default constructor
  ZEventID(const CPid& cpid, int event_id);
  // Constructor for initial event
  ZEventID(bool initial);

  ZEventID(const ZEventID& oth) = default;
  ZEventID(ZEventID&& oth) = default;
  //*** changed
  ZEventID& operator=(const ZEventID& oth) = default;
  ZEventID& operator=(ZEventID&& oth) = default;

  const CPid& cpid() const { return _cpid; }
  int event_id() const { return _event_id; }
  std::size_t hash() const;

  std::string to_string() const;
  void dump() const;

  /* Comparison implements a total order over ZEventIDs. */
  bool operator==(const ZEventID &c) const { return compare(c) == 0; };
  bool operator!=(const ZEventID &c) const { return compare(c) != 0; };
  bool operator<(const ZEventID &c) const { return compare(c) < 0; };
  bool operator<=(const ZEventID &c) const { return compare(c) <= 0; };
  bool operator>(const ZEventID &c) const { return compare(c) > 0; };
  bool operator>=(const ZEventID &c) const { return compare(c) >= 0; };
  int compare(const ZEventID &c) const;

 private:
  /* A complex identifier of the thread */
  CPid _cpid;
  /* Sequential number (within the thread) of this event,
   * The first event of the thread is number 0,
   * the initial event holds event number -1 */
  int _event_id;
  /* Hash */
  std::size_t _hash;

  std::size_t compute_hash() const;
};


namespace std {
  template <>
  struct hash<ZEventID> {
    std::size_t operator()(const ZEventID& k) const {
      return k.hash();
    }
  };
}

#endif // _Z_EVENTID_H_
