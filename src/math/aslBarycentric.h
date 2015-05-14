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


#ifndef ASLBARYCENTRIC_H
#define ASLBARYCENTRIC_H

#include "aslMatrices.h"

namespace asl
{
	/// realize Barycentric algoritms for trianles generated within Kernel 
	/**
		 \ingroup ComplexDataTypes

		 The cordinates does not contain the fist coordinate
	 */
	class Barycentric
	{	
		public:
			std::vector<AVec<>> corners;
			AMatr<> t;
			AMatr<> tInv;
			Barycentric(std::vector<AVec<>> & p);
			Barycentric();
			void init(std::vector<AVec<>> & p);

			AVec<> getCordinates(const AVec<> & p);
			double interpolate(const AVec<> & p, 
			                   const AVec<> & f);
			bool in(const AVec<> & p);
			AVec<> gradient(const AVec<> & f);			
	};
}  //namespace asl

#endif // ASLBARYCENTRIC_H
