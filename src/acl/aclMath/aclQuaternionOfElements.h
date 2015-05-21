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


#ifndef ACLQUATERNIONOFELEMENTS_H
#define ACLQUATERNIONOFELEMENTS_H

#include "aclVectorOfElements.h"

namespace acl
{
	/// The class represents a matrix elements of ::Element
	/**
		 \ingroup ComplexDataTypes
	 */
	class QuaternionOfElements
	{
		private:
			VectorOfElements w; 
			VectorOfElements u; 
		public:	
			QuaternionOfElements();
			
			void setWElement(Element a);
			void setUElement(unsigned int i, Element a);
			const Element getWElement() const;
			const Element getUElement(unsigned int i) const;
			VectorOfElements & getU();
			const VectorOfElements & getU() const;
			VectorOfElements & getW();
			const VectorOfElements & getW() const;
	};

	///function copies the QuaternionOfElements class. 
	/**
		 \relates QuaternionOfElements
	 */
	void copy(const QuaternionOfElements & source, QuaternionOfElements & destination);

	///summ of two matrices
	/**
		 \relates QuaternionOfElements
	 */	
	QuaternionOfElements operator+(QuaternionOfElements & a, QuaternionOfElements & b);
	///difference of two matrices
	/**
		 \relates QuaternionOfElements
	 */	
	QuaternionOfElements operator-(QuaternionOfElements & a,QuaternionOfElements & b);	

	///product of two matrices
	/**
		 \relates QuaternionOfElements
	 */	
	QuaternionOfElements operator*(const QuaternionOfElements & a, const QuaternionOfElements & b);	

	///L2 norm of a quaternion 
	/**
		 \relates QuaternionOfElements
	 */	

	VectorOfElements l2(QuaternionOfElements & a);

	QuaternionOfElements normalize(QuaternionOfElements & a);

	
		
	

}  //namespace acl

#endif // ACLQUATERNIONOFELEMENTS_H
