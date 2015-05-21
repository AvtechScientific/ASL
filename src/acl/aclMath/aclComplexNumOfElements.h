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


#ifndef ACLCOMPLEXNUMOFELEMENTS_H
#define ACLCOMPLEXNUMOFELEMENTS_H

#include "aclVectorOfElements.h"

namespace acl
{
	/// The class represents a matrix elements of ::Element
	/**
		 \ingroup ComplexDataTypes
	 */
	class ComplexNumOfElements
	{
		private:
			VectorOfElements ve; 
		public:	
			ComplexNumOfElements();
			
			void setRe(Element a);
			void setIm(Element a);
			const Element getRe() const;
			const Element getIm() const;
			VectorOfElements & getInternalVector();
			const VectorOfElements & getInternalVector() const;
	};

	///function copies the ComplexNumOfElements class. 
	/**
		 \relates ComplexNumOfElements
	 */
	void copy(const ComplexNumOfElements & source, ComplexNumOfElements & destination);

	///summ of two matrices
	/**
		 \relates ComplexNumOfElements
	 */	
	ComplexNumOfElements operator+(ComplexNumOfElements & a, ComplexNumOfElements & b);
	///difference of two matrices
	/**
		 \relates ComplexNumOfElements
	 */	
	ComplexNumOfElements operator-(ComplexNumOfElements & a,ComplexNumOfElements & b);	

	///product of two matrices
	/**
		 \relates ComplexNumOfElements
	 */	
	ComplexNumOfElements operator*(const ComplexNumOfElements & a, const ComplexNumOfElements & b);	



	
		
	

}  //namespace acl

#endif // ACLCOMPLEXNUMOFELEMENTS_H
