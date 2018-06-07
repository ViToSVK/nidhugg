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
	typedef VCIID AnnotationKeyT;
	typedef int AnnotationValueT; // std::pair<int, CPid>
  typedef std::map<AnnotationKeyT, AnnotationValueT> MappingT;
  typedef MappingT::iterator iterator;
  typedef MappingT::const_iterator const_iterator;
	typedef std::unordered_map<SymAddrSize, std::unordered_set<AnnotationValueT>> ActiveValuesT;

 private:
	
  MappingT mapping;
	ActiveValuesT activevalues;
	ActiveValuesT lockvalues;

 public:
	
  VCAnnotation() = default;
  VCAnnotation(const VCAnnotation&) = default;
  VCAnnotation& operator=(const VCAnnotation&) = default;
  VCAnnotation(VCAnnotation&& a) = default;
  VCAnnotation& operator=(VCAnnotation&& a) = default;

  bool operator==(const VCAnnotation& oth) const {
		return (mapping == oth.mapping) && (lockvalues == oth.lockvalues);
  }

  bool empty() const { return mapping.empty(); }
  size_t size() const { return mapping.size(); }

  iterator begin() { return mapping.begin(); }
  const_iterator begin() const { return mapping.begin(); }
  iterator end() { return mapping.end(); }
  const_iterator end() const { return mapping.end(); }

  void dump() const;
	
  /* *************************** */
  /* MAPPING                     */
  /* *************************** */

 private:

  AnnotationValueT getValue(AnnotationKeyT&& k) const {
    auto it = mapping.find(k);
		assert(it != mapping.end());
    return it->second;
  }

  void add(AnnotationKeyT&& a, const AnnotationValueT& b) {
		auto it = mapping.find(a);
		assert(it == mapping.end());
    mapping.emplace_hint(it, a, b);
  }

  void erase(AnnotationKeyT&& k) {
		assert(defines(k));
    mapping.erase(k);
  }

  bool defines(AnnotationKeyT&& k) const {
    return mapping.find(k) != mapping.end();
  }

 public:

  // We do not keep locks in the mapping
	
  AnnotationValueT getValue(const VCEvent& e) const {
		assert(isRead(a));
    return getValue(AnnotationKeyT(e.cpid, e.instruction_order));
  }

  void add(const VCEvent& a, const AnnotationValueT& b) {
		assert(isRead(a));
    add(AnnotationKeyT(a.cpid, a.instruction_order), b);
	}

	void erase(const VCEvent& ev) {
		assert(isRead(ev));
    mapping.erase(AnnotationKeyT(ev.cpid, ev.instruction_order));
	}

  bool defines(const VCEvent& ev) const {
		assert(isRead(ev));
    return defines(AnnotationKeyT(ev.cpid, ev.instruction_order));
  }	
	

  /* *************************** */
  /* ACTIVE VALUES               */
  /* *************************** */

 private:

	void tryAddActiveValue(const SymAddrSize& ml, const AnnotationValueT& val) {
		auto it = activevalues.find(ml);
    if (it == activevalues.end())
			it = activevalues.emplace_hint(it, ml, std::unordered_set<AnnotationValueT>());
		auto valit = it->second.find(val);
		if (valit == it->second.end())
		  it->second.emplace_hint(valit, val);
	}

	bool isActive(const SymAddrSize& ml, const AnnotationValueT& val) const {
		auto it = activevalues.find(ml);
    if (it == activevalues.end())
			return false;
		auto valit = it->second.find(val);
		return (valit != it->second.end());
	}

 public:

	void tryAddActiveValue(const VCEvent& ev, const AnnotationValueT& val) {
    assert(isRead(ev));
		tryAddActiveValue(ev.ml, val);
	}

	bool isActiveWrite(const VCEvent& ev) const {
    assert(isWrite(ev));
		return isActive(ev.ml, ev.value);
	}
	
  /* *************************** */
  /* LOCK VALUES                 */
  /* *************************** */

 private:
	
	void addLockValue(const SymAddrSize& ml, const AnnotationValueT& val) {
		auto it = lockvalues.find(ml);
    if (it == lockvalues.end())
			it = lockvalues.emplace_hint(it, ml, std::unordered_set<AnnotationValueT>());
		auto valit = it->second.find(val);
		assert(valit == it->second.end());
		it->second.emplace_hint(valit, val);
	}

	void eraseLockValue(const SymAddrSize& ml, const AnnotationValueT& val) {
    auto it = lockvalues.find(ml);
		assert(it != lockvalues.end());
		auto valit = it->second.find(val);
		assert(valit != it->second.end());
		it->second.erase(valit);
	}

	bool definesLockValue(const SymAddrSize& ml, const AnnotationValueT& val) const {
    auto it = lockvalues.find(ml);
		if (it == lockvalues.end())
			return false;
		auto valit = it->second.find(val);
		return (valit != it->second.end());
	}

 public:

	void addLockValue(const VCEvent& ev, const AnnotationValueT& val) {
    assert(isLock(ev));
		addLockValue(ev.ml, val);
	}

	void eraseLockValue(const VCEvent& ev, const AnnotationValueT& val) {
    assert(isLock(ev));
		eraseLockValue(ev.ml, val);
	}

	bool definesLockValue(const VCEvent& ev, const AnnotationValueT& val) const {
    assert(isLock(ev));
		return definesLockValue(ev.ml, val);
	}	

};

#endif // _VC_ANNOTATION_H_
