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


#include "aclMathAlg.h"
#include "aclVectorOfElements.h"
#include <aclGenerators.h>
#include "aclUtilities.h"
#include "acl.h"
#include "DataTypes/aclConstant.h"
#include "../aslUtilities.h"
#include <algorithm>
#include "Kernels/aclKernel.h"

#include <acl/DataTypes/aclArray.h>

#include <acl/aclGenerators.h>

using asl::errorMessage;

using namespace std;

namespace acl
{

	void findNearestDirectionCode(const vector<asl::AVec<>> & directions, 
	                              VectorOfElements v, 
	                              VectorOfElements iDir, 
	                              ExpressionContainer & k)
	{
		acl::TypeID type(getElementType(v));
		VectorOfElements prod(generateVEPrivateVariable(1, type));
		findNearestDirectionCode(directions, v, iDir, prod, k);
	}

	void findNearestDirectionCode(const vector<asl::AVec<>> & directions, 
	                              VectorOfElements v, 
	                              VectorOfElements iDir, 
	                              VectorOfElements scalProd,
	                              ExpressionContainer & k)
	{
		acl::TypeID type(getElementType(v));
		k << (iDir = generateVEConstant(0))
		  << (scalProd = v * generateVEConstant(asl::normalize(directions[0])));
		for (unsigned int i(1); i < directions.size(); ++i){
			VectorOfElements vd(generateVEConstant(asl::normalize(directions[i])));
			k << (iDir = select(iDir, generateVEConstant(i), scalProd < v * vd , type))
			  << (scalProd = max(scalProd, v * vd, type));			
		}
	}

	acl::VectorOfElements differentSign(acl::VectorOfElements d)
	{
		return (acl::fabs(acl::sumOfElements(acl::sign(d))) < (d.size()-.2));
	}

	VectorOfElements generateVEOutOfBoundarySafe(const VectorOfElements & a)
	{
		unsigned int length(getElementsSize(a));
		auto ind(generateVEIndex(length));
		return excerpt(a, (length + ind) % length );
	}

	VectorOfElements generateVEOutOfBoundarySafe(const VectorOfElements & a, 
	                                             const VectorOfElements & outVal)	
	{
		acl::TypeID type(getElementType(a)); 
		unsigned int length(getElementsSize(a));
		auto ind(generateVEIndex(length));
		auto indE(generateVEIndexExt(length));
		return select(excerpt(a, (length + ind) % length), 
		              outVal, 
		              (indE < 0.) || (indE > (length-1)), 
		              type);
	}

	
} // namespace acl
