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

#ifndef _VC_ANNOTATION_H_
#define _VC_ANNOTATION_H_

#include "VCBasis.h"
#include "VCHelpers.h"

#include <memory>
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace std {
  template <>
	struct hash<std::pair<unsigned, unsigned>>
  {
    std::size_t operator()(const std::pair<unsigned, unsigned>& vciid) const
    {
      return
				(hash<unsigned>()(vciid.first) << 12) +
				hash<unsigned>()(vciid.second);
    }
  };
}

class VCAnnotation {
 public:
  typedef std::pair<unsigned, unsigned> VCIID;
	
  enum class Loc {
    LOCAL, REMOTE, ANY
  };

	class Ann {
	public:
		~Ann() { delete goodLocal; }

		Ann(int v, Loc l, std::unordered_set<VCIID>&& gr, bool haslocal, VCIID gl)
			: value(v), loc(l), goodRemote(std::move(gr)),
			goodLocal(haslocal?(new VCIID(gl.first, gl.second)):nullptr) {}

    Ann(const Ann& oth)
			: value(oth.value),
			loc(oth.loc),
			goodRemote(oth.goodRemote),
			goodLocal(oth.goodLocal ?
								new VCIID(oth.goodLocal->first,
													oth.goodLocal->second) :
								nullptr)
				{}
		
	  Ann(Ann&& oth)
			: value(oth.value),
			loc(oth.loc),
			goodRemote(std::move(oth.goodRemote)),
			goodLocal(oth.goodLocal)
				{oth.goodLocal = nullptr; }
		
		const int value;
		const Loc loc;
		const std::unordered_set<VCIID> goodRemote;
		const VCIID * goodLocal;
	};

	using MappingT = std::unordered_map<VCIID, Ann>;
	using LastLockT = std::unordered_map<SymAddrSize, VCIID>;

 private:
	
  MappingT mapping;
	LastLockT lastlock;

 public:
	
  VCAnnotation() = default;
  VCAnnotation(const VCAnnotation&) = default;
  VCAnnotation& operator=(const VCAnnotation&) = default;
  VCAnnotation(VCAnnotation&& a) = default;
  VCAnnotation& operator=(VCAnnotation&& a) = default;

  bool empty() const { return mapping.empty(); }
  size_t size() const { return mapping.size(); }

  typedef MappingT::iterator iterator;
  typedef MappingT::const_iterator const_iterator;	
  iterator begin() { return mapping.begin(); }
  const_iterator begin() const { return mapping.begin(); }
  iterator end() { return mapping.end(); }
  const_iterator end() const { return mapping.end(); }
	
  void dump() const;
	
  /* *************************** */
  /* MAPPING                     */
  /* *************************** */

 public:

	void add(const Node * nd, Ann&& ann) {
		assert(isRead(nd->getEvent()));
		auto key = VCIID(nd->getProcessID(), nd->getEventID());
    auto it = mapping.find(key);
		assert(it == mapping.end());
		assert(ann.loc != Loc::REMOTE || !ann.goodLocal);
		assert(ann.loc != Loc::LOCAL || ann.goodRemote.empty());
		mapping.emplace_hint(it, key, std::move(ann));
	}

	bool defines(const Node * nd) const {
    assert(isRead(nd->getEvent()));
		auto key = VCIID(nd->getProcessID(), nd->getEventID());
		return (mapping.find(key) != mapping.end());
	}

	bool defines(unsigned pid, unsigned evid) const {
		auto key = VCIID(pid, evid);
		return (mapping.find(key) != mapping.end());
	}	

	bool isGood(const Node * writend, const Node * readnd) const {
    assert(isWrite(writend->getEvent()) && isRead(readnd->getEvent()));
		auto key = VCIID(readnd->getProcessID(), readnd->getEventID());
		auto it = mapping.find(key);
		assert(it != mapping.end());
		auto val = VCIID(writend->getProcessID(), writend->getEventID());
		return (it->second.goodRemote.count(val) ||
						(it->second.goodLocal &&
						 *(it->second.goodLocal) == val));
	}

	const Ann& getAnn(const Node * nd) const {
    assert(isRead(nd->getEvent()));
		auto key = VCIID(nd->getProcessID(), nd->getEventID());
		auto it = mapping.find(key);
		assert(it != mapping.end());
		return it->second;
	}

	const Ann& getAnn(unsigned pid, unsigned evid) const {
		auto key = VCIID(pid, evid);
		auto it = mapping.find(key);
		assert(it != mapping.end());
		return it->second;
	}	

  /* *************************** */
  /* LAST LOCK                   */
  /* *************************** */

	void setLastLock(const Node * nd) {
    assert(isLock(nd->getEvent()));
		auto it = lastlock.find(nd->getEvent()->ml);
		if (it == lastlock.end())
			lastlock.emplace_hint(it, nd->getEvent()->ml,
														VCIID(nd->getProcessID(), nd->getEventID()));
		else
			it->second = VCIID(nd->getProcessID(), nd->getEventID());
	}

	std::pair<bool, VCIID> getLastLock(const Node * nd) const {
    assert(isLock(nd->getEvent()));
		auto it = lastlock.find(nd->getEvent()->ml);
		if (it == lastlock.end())
			return {false, VCIID(1337,47)};
		else
			return {true, it->second};
	}

	bool isLastLock(const Node * nd) const {
    assert(isLock(nd->getEvent()));
		auto it = lastlock.find(nd->getEvent()->ml);
		if (it == lastlock.end())
			return false;
		else
			return (it->second != VCIID(nd->getProcessID(), nd->getEventID()));
	}
	
};

namespace std {
  template <>
	struct hash<std::pair<int, VCAnnotation::Loc>>
  {
    std::size_t operator()(const std::pair<int, VCAnnotation::Loc>& val_loc) const
    {
      return hash<int>()(val_loc.first) +
				(size_t) val_loc.second;
    }
  };
}

#endif // _VC_ANNOTATION_H_
