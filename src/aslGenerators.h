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


#ifndef ASLGENERATORS_H
#define ASLGENERATORS_H

#include <memory>
//#include <acl/aclHardware.h>
#include "acl/aclTypes.h"

namespace cl
{
	class CommandQueue;
}

namespace acl
{
	class VectorOfElementsData;
	class VectorOfElements;
	typedef std::shared_ptr<cl::CommandQueue> CommandQueue;	
}

namespace asl
{
	template <typename V> class DataWithGhostNodes;
	typedef DataWithGhostNodes<acl::VectorOfElementsData> DataWithGhostNodesACLData;
	typedef std::shared_ptr<DataWithGhostNodesACLData> SPDataWithGhostNodesACLData;
	typedef DataWithGhostNodes<acl::VectorOfElements> DataWithGhostNodesACL;
	typedef std::shared_ptr<DataWithGhostNodesACL> SPDataWithGhostNodesACL;
	class AbstractDataWithGhostNodes;
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;
	template <typename V> class DataWrapper;
	typedef DataWrapper<acl::VectorOfElementsData> DataWrapperACLData;
	typedef std::shared_ptr<DataWrapperACLData> SPDataWrapperACLData;
	typedef DataWrapper<acl::VectorOfElements> DataWrapperACL;
	typedef std::shared_ptr<DataWrapperACL> SPDataWrapperACL;
	class Block;

	
	/// generates pointer to ACL Data field with \p n components	
	/**
		 \ingroup DataFields
	*/
	template <typename T> SPDataWrapperACLData generateDataContainerACL_SP(const Block &b,
	                                                                       unsigned int n = 1);

	/// generates pointer to ACL Data field with \p n components and \p gN ghost nodes
	/**
		 \ingroup DataFields
	*/
	template <typename T> SPDataWithGhostNodesACLData generateDataContainerACL_SP(const Block &b,
	                                                                              unsigned int n,
	                                                                              unsigned int gN);

	/// generates pointer to ACL Data field with \p n components and \p gN ghost nodes
	/**
		 \ingroup DataFields
	*/
	template <typename T> SPDataWithGhostNodesACLData generateDataContainerACL_SP(const Block &b,
	                                                                              unsigned int n,
	                                                                              unsigned int gN,
	                                                                              acl::CommandQueue queue);
	
	/// generates pointer to ACL Data field with \p n components of type \p t and \p gN ghost nodes
	/**
		 \ingroup DataFields
	*/
	SPDataWithGhostNodesACLData generateDataContainerACL_SP(const Block &b, 
	                                                        acl::TypeID t,
	                                                        unsigned int n, 
	                                                        unsigned int gN, 
	                                                        acl::CommandQueue queue);

	/// generates pointer to ACL Data field with \p n components of type \p t and \p gN ghost nodes
	/**
		 \ingroup DataFields
	*/
	SPDataWithGhostNodesACLData generateDataContainerACL_SP(const Block &b, 
	                                                        acl::TypeID t,
	                                                        unsigned int n, 
	                                                        unsigned int gN);
	
	///	 \ingroup DataFields
	SPDataWrapperACL generateDataContainer_SP(const Block &b, const acl::VectorOfElements & a);	
	
	/// generates pointer to ACL Data field with container \p a and \p gN ghost nodes
	/**
		 \ingroup DataFields
	*/
	SPDataWithGhostNodesACL generateDataContainer_SP(const Block &b, 
	                                                 const acl::VectorOfElements & a, 
	                                                 unsigned int gN);

	
	///	 \ingroup DataFields
	template <typename T> SPDataWithGhostNodesACL generateDataContainerConst_SP(const Block &b, T a, unsigned int gN);	

	/// generates data container without ghost nodes and safe outOfboundary element acces \ingroup DataFields
	SPAbstractDataWithGhostNodes generateDCFullSafe(SPAbstractDataWithGhostNodes d);

	/// generates data container without ghost nodes and safe outOfboundary element acces \ingroup DataFields
	SPAbstractDataWithGhostNodes generateDCFullSafe(SPAbstractDataWithGhostNodes d, double outVal);
	
	

	
}  //namespace acl

#endif // ASLGenerator_H
