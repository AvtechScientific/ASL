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


#include "aslBarycentric.h"
#include "aslMatrices.h"
#include "../aslUtilities.h"

using namespace std;
using asl::errorMessage;

namespace asl
{
	Barycentric::Barycentric(vector<AVec<>> & p)
	{
		init(p);
	}
	
	Barycentric::Barycentric():
		corners(0)
	{
	}

	
	void Barycentric::init(vector<AVec<>> & p)
	{

		if(p.size()-1 != p[0].getSize())
			errorMessage("asl::Barycentric::init: number of points does not corresponds to the dimentionality");
			
		corners=p;

		unsigned int n(p.size()-1);
		t.resize(n, n);
		tInv.resize(n, n);

		for(unsigned int i(1); i < p.size(); ++i)
		{
			t.setRow(i - 1, corners[i]-corners[0]);
		}
		
		tInv = inverseMatrix(t);
	}

	AVec<> Barycentric::getCordinates(const AVec<> & p)
	{
		if(t.getNRow() != p.getSize())
			errorMessage("asl::Barycentric::interpolate: point dimensionality does not corresponds to the triangle dimentionality");
		return tInv * (p - corners[0]);
	}

	double Barycentric::interpolate(const AVec<> & p, 
	                                          const AVec<> & f)
	{
		if(p.getSize()+1 != f.getSize())
			errorMessage("Barycentric::interpolate: number of funciton values does not corresponds to the dimentionality");

		AVec<> fm(subAVec(f, 1, p.getSize()) - AVec<>(p.getSize(), f[0]));
		return f[0] + fm * getCordinates(p); 
	}
		
	bool  Barycentric::in(const AVec<> & p)
	{
		auto x(getCordinates(p));
		auto res(0 <= x[0]);
		for(unsigned int i(1); i < x.getSize(); ++i)
			res = res && ( 0 <= x[i] );
		res = res && (sumOfElements(x) <= 1);
		return res;
	}

	AVec<> Barycentric::gradient(const AVec<> & f)
	{
		if(t.getNRow()+1 != f.getSize())
			errorMessage("Barycentric::gradient: number of funciton values does not corresponds to the dimentionality");
		unsigned int nd(f.getSize()-1);
		AVec<> fm(subAVec(f, 1, nd) - AVec<>(nd,f[0]));
		return fm*tInv; 
	}
		
} // namespace asl
