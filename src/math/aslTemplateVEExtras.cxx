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


#include "aslTemplateVE.h"
#include "acl/acl.h"
#include "acl/aclGenerators.h"
#include <aslUtilities.h>

using acl::VectorOfElements;
using acl::Element;

namespace asl
{

	acl::VectorOfElements cellValues(const TemplateVE & a, unsigned int iEl)
	{
		if(iEl >= a.vto->elementaryCells.size())
			errorMessage("cellValues: iEl larger than the number of cells");
		auto cell(a.vto->elementaryCells[iEl]);
		unsigned int nv(cell.getSize());
		
		acl::VectorOfElements valI(nv);
		for(unsigned int ii(0); ii < nv; ++ii)
			valI[ii] = subVE(a.values,cell[ii])[0];

		return valI;
	}

	
	VectorOfElements gradient(const TemplateVE & a, unsigned int iEl)
	{
		return a.vto->cellGradient(a.values, iEl);
	}
	
	vector<Element> gcGradientAllCells(const TemplateVE & a, 
	                                   vector<VectorOfElements> & values)
	{
		unsigned int nCells(a.vto->elementaryCells.size());
		unsigned int nDim(nD(*a.vectorTemplate));
		auto type(getElementType(a.values));

		values.resize(nCells);
		for(unsigned int i(0); i < nCells; ++i)
			copy(acl::generateVEPrivateVariable(nDim, type),values[i]);

		acl::VectorOfElements code;
		copy(values[0] =  gradient(a, 0), code);
		for(unsigned int i(1); i < nCells; ++i)
			copy(cat(code, values[i] =  gradient(a, i)), 
			     code);			

		return code;
	}
		
		
}// asl
