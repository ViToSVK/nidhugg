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

#include "VCIID.h"
#include "VCHelpers.h"

#include <unordered_map>
#include <unordered_set>

class VCAnnotation {
 public:
  enum class Loc {
    LOCAL, REMOTE, ANY
  };	
	typedef VCIID AnnotationKeyT;
	typedef std::pair<int, Loc> AnnotationValueT;
  typedef std::map<AnnotationKeyT, AnnotationValueT> MappingT;
	typedef std::unordered_map<SymAddrSize, AnnotationKeyT> LastLockT;

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

 private:

  void add(AnnotationKeyT&& a, AnnotationValueT&& b) {
		auto it = mapping.find(a);
		assert(it == mapping.end());
    mapping.emplace_hint(it, a, b);
  }

 public:

  // We DO NOT keep locks in the mapping
  void add(const VCEvent& ev, int val, Loc loc) {
		assert(isRead(ev));
    add(AnnotationKeyT(ev.cpid, ev.instruction_order),
				AnnotationValueT(val, loc));
	}

  /* *************************** */
  /* LAST LOCK                   */
  /* *************************** */

	void setLastLock(const VCEvent& ev) {
    assert(isLock(ev));
		auto it = lastlock.find(ev.ml);
		if (it == lastlock.end())
			lastlock.emplace_hint(it, ev.ml, AnnotationKeyT(ev));
		else
			it->second = AnnotationKeyT(ev);
	}

	std::pair<bool, AnnotationKeyT> getLastLock(const VCEvent& ev) const {
    assert(isLock(ev));
		auto it = lastlock.find(ev.ml);
		if (it == lastlock.end())
			return {false, AnnotationKeyT()};
		else
			return {true, it->second};
	}	
	
};

#endif // _VC_ANNOTATION_H_
