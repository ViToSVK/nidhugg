/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2018 Viktor Toman
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

#ifndef _VCIID_H_
#define _VCIID_H_

#include "VCEvent.h"

class VCIID {
public:
	CPid cpid;
	unsigned instruction_order = 0;

	VCIID() = default;
	VCIID(const CPid& pid, unsigned ord)
		: cpid(pid), instruction_order(ord) {}

	VCIID(const VCIID& oth)
		: cpid(oth.cpid),
			instruction_order(oth.instruction_order) {}

	VCIID& operator=(const VCIID&) = default;

	VCIID(VCIID&& oth)
		: cpid(std::move(oth.cpid)),
			instruction_order(oth.instruction_order) {}

	VCIID& operator=(VCIID&& oth) {
		cpid = std::move(oth.cpid);
		instruction_order = oth.instruction_order;
		return *this;
	}
	
	VCIID(const VCEvent& ev)
		: VCIID(ev.cpid, ev.instruction_order) {}

	// compare the cpid as the last one, because it is the most expensive
	bool operator<(const VCIID& rhs) const {
		return std::tie(instruction_order, cpid)
				< std::tie(rhs.instruction_order, rhs.cpid);
	}

	bool operator==(const VCIID& rhs) const {
		return std::tie(instruction_order, cpid)
				== std::tie(rhs.instruction_order, rhs.cpid);
	}

	bool operator!=(const VCIID& rhs) const {
		return !operator==(rhs);
	}

	void dump() const;
	
};

namespace std {
  template <>
  struct hash<VCIID>
  {
    std::size_t operator()(const VCIID& vciid) const
    {
      return (vciid.instruction_order << 16)
				^ (std::hash<CPid>()(vciid.cpid));
    }
  };
}

#endif // _VCIID_H_
