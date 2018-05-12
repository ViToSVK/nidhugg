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

#endif // _VCIID_H_
