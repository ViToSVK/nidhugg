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
	typedef std::pair<int, CPid> AnnotationValueT;
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

  void intersect(const VCAnnotation& rhs) {
    MappingT new_mapping;
    auto it1 = mapping.begin();
    auto it2 = rhs.mapping.begin();
    auto last1 = mapping.end();
    auto last2 = rhs.mapping.end();

    while (it1 != last1 && it2 != last2) {
			if (*it1 < *it2) {
				++it1;
			} else {
				if (!(*it2 < *it1)) {
						new_mapping.emplace(*it1++);
				}
				++it2;
			}
    }

    new_mapping.swap(mapping);
  }

  size_t size() const { return mapping.size(); }

  const AnnotationValueT *getValue(const AnnotationKeyT& k) const
  {
    auto it = mapping.find(k);
    if (it == mapping.end())
      return nullptr;

    return &it->second;
  }

  const AnnotationValueT *getValue(const VCEvent& e) const
  {
    return getValue(AnnotationKeyT(e.cpid, e.instruction, e.instruction_order));
  }

  // add an annotation (a, b). If there is already an annotation
  // for a, return false.
  bool add(const AnnotationKeyT& a, const AnnotationValueT& b)
  {
    auto it = mapping.find(a);
    if (it != mapping.end())
      return false;

    mapping.emplace_hint(it, a, b);
    return true;
  }

  bool add(AnnotationKeyT&& a, AnnotationValueT&& b)
  {
    auto it = mapping.find(a);
    if (it != mapping.end())
      return false;

    mapping.emplace_hint(it, a, b);
    return true;
  }

  void erase(const AnnotationKeyT& k) {
    mapping.erase(k);
  }

  bool add(const VCEvent& a, const VCEvent& b)
  {
    assert(a.instruction && "Read does not have an instruction");
    return add(AnnotationKeyT(a.cpid, a.instruction, a.instruction_order),
               std::make_pair(b.value, b.cpid));
  }

  bool defines(const AnnotationKeyT& k) const
  {
    return mapping.find(k) != mapping.end();
  }

  bool defines(const VCEvent& ev) const
  {
    return defines(ev.cpid, ev.instruction, ev.instruction_order);
  }

  bool defines(const CPid& c, const llvm::Instruction *i, unsigned ord) const
  {
    return defines(AnnotationKeyT(c, i, ord));
  }

  bool empty() const { return mapping.empty(); }

  iterator begin() { return mapping.begin(); }
  const_iterator begin() const { return mapping.begin(); }
  iterator end() { return mapping.end(); }
  const_iterator end() const { return mapping.end(); }

  void dump() const;
/*
  static PositiveAnnotation getObservationFunction(const std::vector<DCEvent>& trace);
  // returns a positive annotation for an event that forces threads to come
  // to this event
  // \param ev     the event in question
  // \param O      observation function
  // \param trace  trace from which to deduce the past code
  static PositiveAnnotation getPastConeAnnotation(const DCEvent& ev,
                                                  const PositiveAnnotation& O,
                                                  Basis& basis);
*/
private:
  MappingT mapping;
};

#endif // _VC_ANNOTATION_H_
