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
	const llvm::Instruction *instruction = nullptr;
	unsigned instruction_order = 0;

	VCIID() = default;
	VCIID(const CPid& pid, const llvm::Instruction *i, unsigned ord)
		: cpid(pid), instruction(i), instruction_order(ord) {}

	VCIID(const std::tuple<const CPid&, const llvm::Instruction *, unsigned>& tup)
		: cpid(std::get<0>(tup)),
			instruction(std::get<1>(tup)),
			instruction_order(std::get<2>(tup)) {}

	VCIID(const VCIID& oth)
		: cpid(oth.cpid),
			instruction(oth.instruction),
			instruction_order(oth.instruction_order) {}

	VCIID& operator=(const VCIID&) = default;

	VCIID(VCIID&& oth)
		: cpid(std::move(oth.cpid)),
			instruction(oth.instruction),
			instruction_order(oth.instruction_order) {}

	VCIID& operator=(VCIID&& oth) {
		cpid = std::move(oth.cpid);
		instruction = oth.instruction;
		instruction_order = oth.instruction_order;
		return *this;
	}
	
	VCIID(const VCEvent& ev)
		: VCIID(ev.cpid, ev.instruction, ev.instruction_order) {}

	// compare the cpid as the last one, because it is the most expensive
	bool operator<(const VCIID& rhs) const {
		return std::tie(instruction, instruction_order, cpid)
				< std::tie(rhs.instruction, rhs.instruction_order, rhs.cpid);
	}

	bool operator==(const VCIID& rhs) const {
		return std::tie(instruction, instruction_order, cpid)
				== std::tie(rhs.instruction, rhs.instruction_order, rhs.cpid);
	}

	bool operator!=(const VCIID& rhs) const {
		return !operator==(rhs);
	}

	bool operator==(const VCEvent& ev) const {
		return std::tie(instruction, instruction_order, cpid)
				== std::tie(ev.instruction, ev.instruction_order, ev.cpid);
	}

	bool operator!=(const VCEvent& ev) const {
		return !operator==(ev);
	}

	void dump() const;
	
};

namespace std {
  template <>
  struct hash<VCIID>
  {
    std::size_t operator()(const VCIID& vciid) const
    {
      return (vciid.instruction_order >> 8)
				^ (vciid.cpid.hash());
    }
  };
}

#endif // _VCIID_H_
