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


#include "aslDataWrapper.h"
#include "acl/aclUtilities.h"
#include "acl/aclGenerators.h"
#include "acl/aclMath/aclVectorOfElements.h"

namespace asl
{
	template <typename V> bool DataWrapper<V>::checkConsistency() const 
	{
		return compatibleSizes((unsigned int)productOfElements(block.getSize()),
		                        container);
	}

	template bool DataWrapperACL::checkConsistency() const;
	template bool DataWrapperACLData::checkConsistency() const;

	template <class V> const acl::VectorOfElements DataWrapper<V>::getEContainer()
	{
		return acl::VectorOfElements(container);
	}

	template const acl::VectorOfElements DataWrapperACL::getEContainer();
	template const acl::VectorOfElements DataWrapperACLData::getEContainer();
	
	template <> const acl::VectorOfElementsData DataWrapperACL::getDContainer() const
	{
		return acl::VectorOfElementsData();
	}

	template <> const acl::VectorOfElementsData DataWrapperACLData::getDContainer() const
	{
		return container;
	}
	
} // namespace asl
