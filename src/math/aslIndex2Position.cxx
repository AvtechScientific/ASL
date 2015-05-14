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


#include "aslIndex2Position.h"
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include <acl/DataTypes/aclIndexExt.h>
#include <acl/DataTypes/aclConstant.h>

namespace asl
{	
	Index2PositionDiscreteACL::Index2PositionDiscreteACL(const Block & b, bool putSize):
		initPosition(b.getSize().getSize()),
		position(acl::generateVEPrivateVariable<int>(b.getSize().getSize())),
		positionWithInit(b.getSize().getSize())		
	{
		unsigned int nD(b.getSize().getSize());
		unsigned int size(putSize ? productOfElements (b.getSize()) : 0);
		if (nD==2){
			using namespace acl::elementOperators;
			using acl::elementOperators::operator/;
			positionWithInit[0] =acl::Element(new acl::IndexExt(size))/
								acl::Element(new acl::Constant<int>(b.getSize()[1]));
			positionWithInit[1] =acl::Element(new acl::IndexExt(size))%
								acl::Element(new acl::Constant<int>(b.getSize()[1]));
			copy(position=positionWithInit,initPosition);
		}
		if (nD==3){
			using namespace acl::elementOperators;
			int s12(b.getSize()[1]*b.getSize()[2]);
			using acl::elementOperators::operator/;
			positionWithInit[0] =acl::Element(new acl::IndexExt(size))/
								acl::Element(new acl::Constant<int>(s12));
			positionWithInit[1] =acl::Element(new acl::IndexExt(size))%
								acl::Element(new acl::Constant<int>(s12))/
								acl::Element(new acl::Constant<int>(b.getSize()[2]));
			positionWithInit[2] =acl::Element(new acl::IndexExt(size))%
								acl::Element(new acl::Constant<int>(s12))%
								acl::Element(new acl::Constant<int>(b.getSize()[2]));
			copy(position=positionWithInit,initPosition);
		}
	}

	Index2PositionACL::Index2PositionACL(const Block & b, 
	                                     acl::TypeID type, 
	                                     bool putSize):
		initPosition(b.getSize().getSize()),
		position(acl::generateVEPrivateVariable(b.getSize().getSize(),type)),
		positionWithInit(b.getSize().getSize())		
	{
		unsigned int nD(b.getSize().getSize());
		unsigned int size(putSize ? productOfElements (b.getSize()) : 0);
		if (nD==2){
			using namespace acl::elementOperators;
			using acl::elementOperators::operator/;
			positionWithInit[0] =acl::Element(new acl::IndexExt(size))/
								acl::Element(new acl::Constant<int>(b.getSize()[1]));
			positionWithInit[1] =acl::Element(new acl::IndexExt(size))%
								acl::Element(new acl::Constant<int>(b.getSize()[1]));
			copy(acl::convert(type, 
			                  acl::generateVEConstant(b.position) + 
			                  acl::generateVEConstant(b.dx) * 
			                  acl::convert(type,positionWithInit),
			                  false), 
			     positionWithInit);
			copy(position=positionWithInit, initPosition);
		}
		if (nD==3){
			using namespace acl::elementOperators;
			int s12(b.getSize()[1]*b.getSize()[2]);
			using acl::elementOperators::operator/;
			positionWithInit[0] =acl::Element(new acl::IndexExt(size))/
								acl::Element(new acl::Constant<int>(s12));
			positionWithInit[1] =acl::Element(new acl::IndexExt(size))%
								acl::Element(new acl::Constant<int>(s12))/
								acl::Element(new acl::Constant<int>(b.getSize()[2]));
			positionWithInit[2] =acl::Element(new acl::IndexExt(size))%
								acl::Element(new acl::Constant<int>(s12))%
								acl::Element(new acl::Constant<int>(b.getSize()[2]));
			copy(acl::convert(type, 
   			                  acl::generateVEConstant(b.position) + 
			                  acl::generateVEConstant(b.dx) * 
			                  acl::convert(type,positionWithInit), 
			                  false),
			     positionWithInit);			
			copy(position=positionWithInit, initPosition);
		}
	}

		
}// asl
