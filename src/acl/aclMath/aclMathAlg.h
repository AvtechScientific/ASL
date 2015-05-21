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


#ifndef ACLMATHALG_H
#define ACLMATHALG_H

#include "aclVectorOfElementsDef.h"
#include "math/aslVectorsDynamicLength.h"

using namespace std;

namespace acl
{
	class ExpressionContainer;
	
	/// generates code for finding nearest direction from given directions set \p directions
	/**
		 \param directions vectors set defining directions
		 \param v the vector for which directions is seeked
		 \param iDir the number of vector is stored here
		 \param k the kernel for the generated code
	*/
	void findNearestDirectionCode(const vector<asl::AVec<>> & directions, 
	                           VectorOfElements v, 
	                           VectorOfElements iDir, 
	                           ExpressionContainer & k);

	/// generates code for finding nearest direction from given directions set \p directions
	/**
		 \param directions vector set defining directions the directions are normalized internaly
		 \param v the vector for which directions is seeked
		 \param iDir the number of vector is stored here
		 \param scalProd the scalar product of the vector with the corresponding normalized direction vector is stored here 
		 \param k the kernel for the generated code
	*/
	void findNearestDirectionCode(const vector<asl::AVec<>> & directions, 
	                           VectorOfElements v, 
	                           VectorOfElements iDir,
	                           VectorOfElements scalProduct,
	                           ExpressionContainer & k);

	/// generate expresion returning true if elements of v have different signs or one of them zerow 
	VectorOfElements differentSign(VectorOfElements v);

	/// generates Vector of elements wraping the \p a in order to avoid out of boundary acces
	VectorOfElements generateVEOutOfBoundarySafe(const VectorOfElements & a);	

	/// generates Vector of elements wraping the \p a in order to avoid out of boundary acces with given out of boundary value
	VectorOfElements generateVEOutOfBoundarySafe(const VectorOfElements & a, 
	                                             const VectorOfElements & outVal);	

	
}  //namespace acl

#endif // ACLMATHALG_H
