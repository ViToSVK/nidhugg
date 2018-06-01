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

class VCAnnotation {
 public:
	typedef VCIID AnnotationKeyT;
	typedef int AnnotationValueT; // std::pair<int, CPid>
  typedef std::map<AnnotationKeyT, AnnotationValueT> MappingT;
  typedef MappingT::iterator iterator;
  typedef MappingT::const_iterator const_iterator;

  const MappingT& getMapping() const { return mapping; }

  VCAnnotation() = default;
  VCAnnotation(const VCAnnotation&) = default;
  VCAnnotation& operator=(const VCAnnotation&) = default;
  VCAnnotation(VCAnnotation&& a) = default;
  VCAnnotation& operator=(VCAnnotation&& a) = default;

  bool operator==(const VCAnnotation& oth) const {
      return mapping == oth.mapping;
  }

  size_t size() const { return mapping.size(); }

  const AnnotationValueT *getValue(const AnnotationKeyT& k) const
  {
    auto it = mapping.find(k);
		assert(it != mapping.end());
    return &(it->second);
  }

  const AnnotationValueT *getValue(const VCEvent& e) const
  {
    return getValue(AnnotationKeyT(e.cpid, e.instruction_order));
  }

  void add(const AnnotationKeyT& a, const AnnotationValueT& b)
  {
		assert(mapping.find(a) == mapping.end());
    mapping.emplace(a, b);
  }

  void add(AnnotationKeyT&& a, AnnotationValueT&& b)
  {
		assert(mapping.find(a) == mapping.end());
		mapping.emplace(a, b);
  }

  void add(const VCEvent& a, const VCEvent& b)
  {
		assert(isRead(a) && isWrite(b));
    add(AnnotationKeyT(a.cpid, a.instruction_order),
				b.value); //  std::make_pair(b.value, b.cpid)
  }
	
  void erase(const AnnotationKeyT& k) {
		assert(defines(k));
    mapping.erase(k);
  }

	void erase(const VCEvent& ev) {
    mapping.erase(AnnotationKeyT(ev.cpid, ev.instruction_order));
	}

  bool defines(const AnnotationKeyT& k) const
  {
    return mapping.find(k) != mapping.end();
  }

  bool defines(const VCEvent& ev) const
  {
    return defines(AnnotationKeyT(ev.cpid, ev.instruction_order));
  }

  bool empty() const { return mapping.empty(); }

  iterator begin() { return mapping.begin(); }
  const_iterator begin() const { return mapping.begin(); }
  iterator end() { return mapping.end(); }
  const_iterator end() const { return mapping.end(); }

  void dump() const;

private:
  MappingT mapping;
};

#endif // _VC_ANNOTATION_H_
