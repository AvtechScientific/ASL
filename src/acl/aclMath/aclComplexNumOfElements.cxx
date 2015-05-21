/*
 * Advanced Simulation Library <http://asl.org.il>
 * 
 * Copyright 2015 Avtech Scientific <http://avtechscientific.com>
 *
 *
 * This file is part of Advanced Simulation Library (ASL).
 *
 * ASL is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, version 3 of the License.
 *
 * ASL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with ASL. If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "aclComplexNumOfElements.h"
#include "aclUtilities.h"
#include "acl.h"
#include "../aslUtilities.h"
#include <algorithm>

using namespace std;
using asl::errorMessage;

namespace acl
{
	ComplexNumOfElements::ComplexNumOfElements():
		ve(2)
	{
	}

	
	void ComplexNumOfElements::setRe(Element a)
	{
		ve[0]=a;
	}

	void ComplexNumOfElements::setIm(Element a)
	{
		ve[1]=a;
	}
		
	const Element ComplexNumOfElements::getRe() const
	{
		return ve[0];
	}

	const Element ComplexNumOfElements::getIm() const
	{
		return ve[1];
	}
	
	VectorOfElements & ComplexNumOfElements::getInternalVector()
	{
		return ve;
	}

	const VectorOfElements & ComplexNumOfElements::getInternalVector() const
	{
		return ve;
	}

		
	void copy(const ComplexNumOfElements & source,ComplexNumOfElements  destination)
	{
		copy(source.getInternalVector(),destination.getInternalVector());		
	}

	ComplexNumOfElements operator+(ComplexNumOfElements & a,ComplexNumOfElements & b)
	{
		ComplexNumOfElements m;
		copy(a.getInternalVector()+b.getInternalVector(),m.getInternalVector());
		return m;			
	}

	ComplexNumOfElements operator-(ComplexNumOfElements & a,ComplexNumOfElements & b)
	{
		ComplexNumOfElements m;
		copy(a.getInternalVector()-b.getInternalVector(),m.getInternalVector());
		return m;
	}
		
	ComplexNumOfElements operator*(const ComplexNumOfElements & a, const ComplexNumOfElements & b)
	{
		ComplexNumOfElements m;
		using namespace elementOperators;
		m.setRe(a.getRe()*b.getRe()-a.getIm()*b.getIm());
		m.setIm(a.getRe()*b.getIm()+a.getIm()*b.getRe());
		return m;			
	}

} // namespace acl
