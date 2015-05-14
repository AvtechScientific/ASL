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


#include "aclQuaternionOfElements.h"
#include "aclVectorOfElementsOperations.h"
#include "aclGenerators.h"
#include "aclElementBase.h"

using namespace std;
using asl::errorMessage;

namespace acl
{
	QuaternionOfElements::QuaternionOfElements():
		w(1),
		u(3)
	{
	}

	
	void QuaternionOfElements::setWElement(Element a)
	{
		w[0]=a;
	}

	void QuaternionOfElements::setUElement(unsigned int i, Element a)
	{
		u.at(i)=a;
	}
		
	const Element QuaternionOfElements::getWElement() const
	{
		return w[0];
	}

	const Element QuaternionOfElements::getUElement(unsigned int i) const
	{
		return u.at(i);
	}
	
	VectorOfElements & QuaternionOfElements::getW()
	{
		return w;
	}

	VectorOfElements & QuaternionOfElements::getU()
	{
		return u;
	}

	const VectorOfElements & QuaternionOfElements::getW() const
	{
		return w;
	}

	const VectorOfElements & QuaternionOfElements::getU() const
	{
		return u;
	}
		
		
	void copy(const QuaternionOfElements & source,QuaternionOfElements  destination)
	{
		copy(source.getW(),destination.getW());		
		copy(source.getU(),destination.getU());		
	}

	QuaternionOfElements operator+(QuaternionOfElements & a,QuaternionOfElements & b)
	{
		QuaternionOfElements m;
		copy(a.getW()+b.getW(),m.getW());
		copy(a.getU()+b.getU(),m.getU());
		return m;			
	}

	QuaternionOfElements operator-(QuaternionOfElements & a,QuaternionOfElements & b)
	{
		QuaternionOfElements m;
		copy(a.getW()-b.getW(),m.getW());
		copy(a.getU()-b.getU(),m.getU());
		return m;
	}
		
	QuaternionOfElements operator*(const QuaternionOfElements & a, const QuaternionOfElements & b)
	{
		QuaternionOfElements m;
		copy(a.getW()*b.getW()-a.getU()*b.getU(),m.getW());
		copy(a.getW() * b.getU() + b.getW()*a.getU() +
		     crossProduct(a.getU(), b.getU()), m.getU());
		return m;			
	}

	VectorOfElements l2(QuaternionOfElements & a)
	{
		return l2(a.getW())+l2(a.getU());
	}
		
	QuaternionOfElements normalize(QuaternionOfElements & a)
	{
		VectorOfElements il(generateVEPrivateVariable(1u,a.getW()[0]->getTypeID()));
		VectorOfElements ileq(1);
		copy(il=1./sqrt(l2(a)),ileq);
		QuaternionOfElements m;
		copy(a.getW()*ileq,m.getW());
		copy(a.getU()*il,m.getU());
		return m;
	}	

} // namespace acl
