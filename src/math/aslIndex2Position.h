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


#ifndef INDEX2POSITION_H_INCLUDED
#define INDEX2POSITION_H_INCLUDED

#include <data/aslBlocks.h>
#include <acl/aclMath/aclVectorOfElements.h>


namespace asl
{

	class Index2PositionACL
	{
		public:
			Index2PositionACL(const Block &b, 
			                  acl::TypeID type=acl::TYPE_DOUBLE, 
			                  bool putSize=true);
			acl::VectorOfElements initPosition;
			acl::VectorOfElements position;
			acl::VectorOfElements positionWithInit;
	};

	class Index2PositionDiscreteACL
	{
		public:
			/// \p putSize is a flag defining the size of the generated elements in case of false value the size is 0 
			Index2PositionDiscreteACL(const Block &b, bool putSize=true);
			acl::VectorOfElements initPosition;
			acl::VectorOfElements position;		
			acl::VectorOfElements positionWithInit;
	};


}// asl

#endif // INDEX2POSITION_H_INCLUDED
