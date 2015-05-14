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


#ifndef ACLBARYCENTRIC_H
#define ACLBARYCENTRIC_H

#include "aclMatrixOfElements.h"

namespace acl
{
	/// realizes Barycentric algorithms for triangles generated within Kernel 
	/**
		 \ingroup ComplexDataTypes

		 The coordinates do not contain the first coordinate
	 */
	class Barycentric
	{	
		public:
			vector<acl::VectorOfElements> corners;
			MatrixOfElements t;
			MatrixOfElements tInv;
			VectorOfElements initTInv;
			Barycentric(vector<acl::VectorOfElements> & p);
			Barycentric();
			void init(vector<VectorOfElements> & p);

			VectorOfElements getCordinates(const VectorOfElements & p);
			VectorOfElements interpolate(const VectorOfElements & p, 
			                             const VectorOfElements & f);
			VectorOfElements in(const VectorOfElements & p);
			VectorOfElements gradient(const VectorOfElements & f);
			
	};
}  //namespace acl

#endif // ACLBARYCENTRIC_H
