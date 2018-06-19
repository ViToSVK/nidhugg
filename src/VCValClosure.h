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

#ifndef _VC_VALCLOSURE_H
#define _VC_VALCLOSURE_H

#include "VCGraphVclock.h"

class VCValClosure {
 public:

	typedef std::pair<int, VCAnnotation::Loc> AnnotationValueT;
	
	VCValClosure(const VCGraphVclock& gr, const std::unordered_map
						 <const Node *, AnnotationValueT>& vf)
	: graph(gr), valFunction(vf) {}

	VCValClosure(const VCValClosure& oth) = default;
	VCValClosure& operator=(VCValClosure& oth) = delete;
	VCValClosure(VCValClosure&& oth) = default;
	VCValClosure& operator=(VCValClosure&& oth) = delete;

	//

  void prepareOne(const PartialOrder& po, const Node * readnd);
	
	void prepare(const PartialOrder& po, const Node * newread);  	
	
	void valClose(const PartialOrder& po, const Node * newread,
								const std::pair<int, VCAnnotation::Loc> * newval);

	//

	const VCGraphVclock& graph;
	
	const std::unordered_map
	<const Node *, AnnotationValueT>& valFunction;
	
	bool result;
	
	std::unordered_map
	<SymAddrSize, std::vector<const Node *>> wNonroot;
	
	std::unordered_map
	<const Node *, const std::vector<const Node *>&> wRem;

	std::unordered_map
	<const Node *, std::pair<int, int>> wBounds;

	std::unordered_map
	<const Node *, const Node *> wLoc;
};
	
#endif // _VC_VALCLOSURE_H
