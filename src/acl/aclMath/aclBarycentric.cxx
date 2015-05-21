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


#include "aclBarycentric.h"
#include "aclMatrixOfElements.h"
#include "aclVectorOfElements.h"
#include "aclGenerators.h"
#include "../aslUtilities.h"

using namespace std;
using asl::errorMessage;

namespace acl
{
	Barycentric::Barycentric(vector<acl::VectorOfElements> & p)
	{
		init(p);
	}
	
	Barycentric::Barycentric():
		corners(0)
	{
	}
	
	void Barycentric::init(vector<VectorOfElements> & p)
	{
		acl::TypeID type(getElementType(p[0]));

		if (p.size()-1 != p[0].size())
			errorMessage("Barycentric::init: number of points does not corresponds to the dimentionality");
			
		corners.resize(p.size());
		for (unsigned int i(0); i < p.size(); ++i)
			copy(p[i],corners[i]);

		unsigned int n(p.size()-1);
		t.resize(n, n);
		if (n != tInv.getNRows())
		{
			tInv.resize(n, n);
			copy(generateVEPrivateVariable(n*n, type), tInv.getInternalVector());
		}

		for (unsigned int i(1); i < p.size(); ++i)
		{
			t.setRow(i - 1, corners[i]-corners[0]);
		}
		
		copy(gcMatrixInversion(t,tInv), initTInv);
	}

	VectorOfElements Barycentric::getCordinates(const VectorOfElements & p)
	{
		if(t.getNRows() != p.size())
			errorMessage("Barycentric::interpolate: point dimensionality does not corresponds to the triangle dimentionality");
		return tInv * (p - corners[0]);
	}

	VectorOfElements Barycentric::interpolate(const VectorOfElements & p, 
	                                          const VectorOfElements & f)
	{
		if(p.size()+1 != f.size())
			errorMessage("Barycentric::interpolate: number of funciton values does not corresponds to the dimentionality");

		VectorOfElements fm(subVE(f, 1, p.size()) - catN(subVE(f, 0), p.size()));
		return subVE(f,0) + fm * getCordinates(p); 
	}
		
	VectorOfElements  Barycentric::in(const VectorOfElements & p)
	{
		auto x(getCordinates(p));
		auto zero(generateVEConstant(0));
		auto res(zero <= subVE(x,0));
		for(unsigned int i(1); i < x.size(); ++i)
			copy(res && ( zero <= subVE(x,i) ), res);
		copy(res && ( sumOfElements(x) <= generateVEConstant(1)), res);
		return res;
	}

	VectorOfElements Barycentric::gradient(const VectorOfElements & f)
	{
		if(t.getNRows()+1 != f.size())
			errorMessage("Barycentric::gradient: number of funciton values does not corresponds to the dimentionality");
		unsigned int nd(f.size()-1);
		VectorOfElements fm(subVE(f, 1, nd) - catN(subVE(f, 0), nd));
		return fm*tInv; 
	}
		
} // namespace acl
