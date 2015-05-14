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


#include "aslTemplatesExtras.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include "acl/aclGenerators.h"
#include "aslBarycentric.h"

namespace asl
{
	
	VTObjects::VTObjects (const VectorTemplate* vt, 
	                      const std::vector<unsigned int> & ep1, 
	                      const std::vector<unsigned int> & ep2,
	                      const std::vector<AVec<int>> & elCells):
		vt(vt),
		edgePoint1(ep1),
		edgePoint2(ep2),
		elementaryCells(elCells)
	{
		initCellMatrices();
	}

	void VTObjects::getCellPoints(unsigned int ic, vector<AVec<>> & points) const
	{
		points.resize(elementaryCells[ic].getSize());
		for(unsigned int i(0); i < points.size(); ++i)
		{
			points[i] = vt->vectors[elementaryCells[ic][i]];
		}
	}
		
	void VTObjects::initCellMatrices()
	{
		unsigned int n(elementaryCells.size());
		cellMatrices.resize(n);		
		for(unsigned int i(0); i < n; ++i)
		{
			vector<AVec<>> points;
			getCellPoints(i, points);
			asl::Barycentric b(points);
			cellMatrices[i] = b.tInv;
		}
	}

	acl::VectorOfElements VTObjects::cellGradient(const acl::VectorOfElements & val, 
	                                              unsigned int i) const
	{
		if(i >= elementaryCells.size())
			errorMessage("VTObjects::cellGradient: i larger than the number of cells");
		auto cell(elementaryCells[i]);
		unsigned int nv(cell.getSize());
		
		acl::VectorOfElements val0(subVE(val,cell[0]));
		acl::VectorOfElements valI(nv-1);
		for(unsigned int ii(1); ii < nv; ++ii)
			valI[ii-1] = (subVE(val,cell[ii])-val0)[0];
		
		return valI * acl::generateMEConstant(cellMatrices[i]);
	}

	double edgeWeight(const VTObjects & vto, unsigned int iEl, unsigned int i, unsigned int j)
	{
		auto & cell(vto.elementaryCells[iEl]);
		auto v(vto.vt->vectors[cell[i]]-vto.vt->vectors[cell[j]]);
		unsigned int nDir(vto.vt->vectors.size());
		unsigned int ii(0);
		bool test(true);
		for(; test && ii < nDir; ++ii)
		{
			test = (l2(v) != l2(vto.vt->vectors[ii]));
		}
		double res(0);
		if (!test)
			res=vto.vt->laplasCoefs[ii];
		return res;
	}
		
	const VTObjects & d2q5Objs()
	{
		static VTObjects vto(&d2q5(), 
	                   {0,0,0,0,1,1,3,3},
	                   {1,2,3,4,2,4,2,4},
	                   {makeAVec(0,1,2), makeAVec(0,2,3), makeAVec(0,3,4), makeAVec(0,4,1)});
		return vto;
	}

	const VTObjects & d2q9Objs()
	{
		static VTObjects vto(&d2q9(), 
	                   {0,0,0,0,0,0,0,0,1,5,2,6,3,7,4,8},
	                   {1,2,3,4,5,6,7,8,5,2,6,3,7,4,8,1},
	                   {makeAVec(1,0,8), makeAVec(1,0,5), 
						makeAVec(2,0,5), makeAVec(2,0,6),
						makeAVec(3,0,6), makeAVec(3,0,7),
						makeAVec(4,0,7), makeAVec(4,0,8)});
		return vto;
	}

		
	const VTObjects & d3q7Objs()
	{
		static VTObjects vto(&d3q7(), 
	                {0,0,0,0,0,0,1,1,1,1,4,4,4,4},
	                {1,2,3,4,5,6,2,3,5,6,2,3,5,6},
	                {makeAVec(0,1,2,3), makeAVec(0,2,4,3), makeAVec(0,4,5,3),
					 makeAVec(0,5,1,3), makeAVec(0,2,1,6), makeAVec(0,4,2,6),
					 makeAVec(0,5,4,6), makeAVec(0,1,5,6)});
		return vto;
	}

	const VTObjects & d3q15Objs()
	{
		static VTObjects vto(&d3q15(), 
	                {0,0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0,1,1,1, 1,2,2, 2, 2,3,3, 3, 3, 4, 4, 4, 4,5, 5, 5, 5,6, 6, 6, 6},
	                {1,2,3,4,5,6,7,8,9,10,11,12,13,14,7,8,9,10,7,8,11,12,7,9,11,13,11,12,13,14,9,10,13,14,8,10,12,14},
	                {makeAVec(0,1, 7, 9), makeAVec(0,1, 9,10), makeAVec(0,1,10, 8), makeAVec(0,1, 8, 7),
					 makeAVec(0,2, 7, 8), makeAVec(0,2, 8,12), makeAVec(0,2,12,11), makeAVec(0,2,11, 7),
					 makeAVec(0,3, 7,11), makeAVec(0,3,11,13), makeAVec(0,3,13, 9), makeAVec(0,3, 9, 7),
					 makeAVec(0,4,11,13), makeAVec(0,4,13,14), makeAVec(0,4,14,12), makeAVec(0,4,12,11),
					 makeAVec(0,5, 9,10), makeAVec(0,5,10,14), makeAVec(0,5,14,13), makeAVec(0,5,13, 9),
					 makeAVec(0,6, 8,12), makeAVec(0,6,12,14), makeAVec(0,6,14,10), makeAVec(0,6,10, 8)});
		return vto;
	}

	const VTObjects & d3q19Objs()
	{
		static VTObjects vto(&d3q19(), 
	                {0,0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1, 1,1, 1, 2, 2,2, 2, 3, 3, 3, 3,4, 4, 4, 4,5, 5,5, 5, 6, 6, 6, 6, 7,15, 8,16,10,11, 7,12,11,18,14,15,10,18, 9,17, 9,14, 8,13,12,17,13,16},
	                {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,7,15,8,16,10,11,7,12,11,18,14,15,9,17,10,18,9,14,8,13,12,17,13,16,15, 8,16, 7,11, 7,12,10,18,14,15,11,18, 9,17,10,14, 8,13, 9,17,13,16,12},
	                {makeAVec(0,1, 7,15), makeAVec(0,1,15, 8), makeAVec(0,1, 8,16), makeAVec(0,1,16, 7),
					 makeAVec(0,2,10,11), makeAVec(0,2,11, 7), makeAVec(0,2, 7,12), makeAVec(0,2,12,10),
					 makeAVec(0,3,11,18), makeAVec(0,3,18,14), makeAVec(0,3,14,15), makeAVec(0,3,15,11),
					 makeAVec(0,4,10,18), makeAVec(0,4,18, 9), makeAVec(0,4, 9,17), makeAVec(0,4,17,10),
					 makeAVec(0,5, 9,14), makeAVec(0,5,14, 8), makeAVec(0,5, 8,13), makeAVec(0,5,13, 9),
					 makeAVec(0,6,12,17), makeAVec(0,6,17,13), makeAVec(0,6,13,16), makeAVec(0,6,16,12),
					 makeAVec(0,7,16,12), makeAVec(0,7,11,15), makeAVec(0,8,16,13), makeAVec(0,8,15,14),
					 makeAVec(0,9,14,18), makeAVec(0,9,13,17), makeAVec(0,10,18,11), makeAVec(0,10,17,12)});
		return vto;
	}
		
	const VTObjects* vtObject(const VectorTemplate *vt)
	{
		const VTObjects* vto(NULL);

		if (vt == &d2q5()) 
				vto=&d2q5Objs();
		if (vt == &d2q9()) 
				vto=&d2q9Objs();
		if (vt == &d3q7())
				vto=&d3q7Objs();
		if (vt == &d3q15())
				vto=&d3q15Objs();
		if (vt == &d3q19())
				vto=&d3q19Objs();

//		if	(vto == NULL)
//				errorMessage("vtObject:  VTObjects* is undefined for the given VectorTemplate");
		return vto;
	}
		
}// asl
