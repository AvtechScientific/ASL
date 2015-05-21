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


#include "aslDataWithGhostNodes.h"
#include "acl/aclGenerators.h"
#include "acl/aclMath/aclVectorOfElements.h"

namespace asl
{

	template<class V> acl::VectorOfElements DataWithGhostNodes<V>::getSubContainer ()
	{
		return acl::generateVESubElements(dw.getEContainer(),getSubContainerSize(),getSubContainerOffset());
	}

	template acl::VectorOfElements DataWithGhostNodesACLData::getSubContainer ();
	template acl::VectorOfElements DataWithGhostNodesACL::getSubContainer ();

	template<class V> const acl::VectorOfElements DataWithGhostNodes<V>::getEContainer ()
	{
		return dw.getEContainer();
	}

	template const acl::VectorOfElements DataWithGhostNodesACLData::getEContainer ();
	template const acl::VectorOfElements DataWithGhostNodesACL::getEContainer ();

	template<class V>const  acl::VectorOfElementsData DataWithGhostNodes<V>::getDContainer () const
	{
		return dw.getDContainer();
	}

	template const acl::VectorOfElementsData DataWithGhostNodesACLData::getDContainer () const;
	template const acl::VectorOfElementsData DataWithGhostNodesACL::getDContainer () const;

	SPDataWithGhostNodesACLData clone(SPDataWithGhostNodesACLData d)
	{
		auto newd(std::make_shared<DataWithGhostNodesACLData>(d->getInternalBlock(),
		                                                      d->getGhostBorder()));
		newd->setContainer(clone(d->getDContainer()));
		return newd;
	}

	SPDataWithGhostNodesACLData clone(SPDataWithGhostNodesACLData d, unsigned int n)
	{
		auto newd(std::make_shared<DataWithGhostNodesACLData>(d->getInternalBlock(),
		                                                      d->getGhostBorder()));
		newd->setContainer(clone(d->getDContainer(), n));
		return newd;
	}
	
}// asl
