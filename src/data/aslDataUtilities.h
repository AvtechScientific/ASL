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


#ifndef ASLDATAUTILITIES_H
#define ASLDATAUTILITIES_H

#include "aslDataWrapper.h"
#include "acl/aclMath/aclVectorOfElementsDef.h"

namespace acl
{
	class Kernel;
}
namespace asl
{

	/// Uploads \p source from global to local memory in portion
	/// defined by \p size with group size \p groupSize.
	/// Returns destination (local VectorOfElements).
	/// \ingroup DataUtilities
	acl::VectorOfElements uploadToLocalMem(AbstractData & source,
	                                       const AVec<int> & size,
	                                       unsigned int groupSize,
	                                       acl::Kernel & kernel);

	/// generates DataWraper with points placed inside a widow
	/**
		 \param a the point corresponding to 0 coner
		 \param b the point corresponding to diagonal to 0 one coner
		 \ingroup DataUtilities

		 \todo errorMessages
	*/
	SPDataWrapperACL generateSubData(SPDataWrapperACL d, AVec<int> a, AVec<int> b);
	
}
#endif // ASLDATAUTILITIES_H

