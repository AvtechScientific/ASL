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


#ifndef ASLDATAWRAPPER_H
#define ASLDATAWRAPPER_H

#include "aslBlocks.h"
#include "../acl/aclMath/aclVectorOfElementsDef.h"
#include <iostream>
#include <fstream>

/**
 \defgroup DataFields Data Fileds
 \ingroup Numerics
 */

namespace acl
{
	void copy(const vector<Element> & source,
	          VectorOfElements & destination);
	void copy(const VectorOfElementsData & source,
	          VectorOfElementsData & destination);
}


using  namespace std;

namespace asl
{

	/// \ingroup DataFields
	class AbstractData
	{
		protected:
			Block block;
			inline AbstractData();
			inline explicit AbstractData(const Block & b);
		public:	
			virtual const acl::VectorOfElements getEContainer() = 0;	
			virtual const acl::VectorOfElementsData getDContainer() const = 0;	
			inline const Block & getBlock() const;
			inline void setBlock(const Block & b);
	};

	typedef shared_ptr<AbstractData> SPAbstractData;

		
    /// Class defines a folding rule into a 1D/2D/3D bulk 
    /** 
		 \param V is type of the container
		 \ingroup DataFields

		 \todo resolve consistency issue with setBlock()
		 
    */
	template <typename V> class DataWrapper: public AbstractData
	{
		protected:
			V container;
			virtual bool checkConsistency() const;
		public:
			inline DataWrapper();
			inline explicit DataWrapper(const Block & b);
			inline DataWrapper(DataWrapper & d);
			virtual const acl::VectorOfElements getEContainer();	
			virtual const acl::VectorOfElementsData getDContainer() const;	
			inline V & getContainer();	
			inline void setContainer(const V & cont);
	};

	typedef DataWrapper<acl::VectorOfElements> DataWrapperACL;
	typedef DataWrapper<acl::VectorOfElementsData> DataWrapperACLData;

	typedef shared_ptr<DataWrapperACL> SPDataWrapperACL;
	typedef shared_ptr<DataWrapperACLData> SPDataWrapperACLData;

	
// ---------------------------- Implementation ---------------------------

	AbstractData::AbstractData()
	{
	}
	
	AbstractData::AbstractData(const Block & b): 
		block(b)
	{
	}

	template <typename V> DataWrapper<V>::DataWrapper():
		AbstractData()
	{
	}
	
	template <typename V> DataWrapper<V>::DataWrapper(const Block & b): 
		AbstractData(b)
	{		
	}

	template <typename V> DataWrapper<V>::DataWrapper(DataWrapper & d): 
		AbstractData(d.block),
		container(d.container)
	{		
	}
		
	template <typename V> inline V & DataWrapper<V>::getContainer()
	{
		return container;
	}
	
	template <typename V> inline  void DataWrapper<V>::setContainer(const V & cont)
	{
		copy(cont, container);
	}

	inline const Block & AbstractData::getBlock() const
	{
		return block;
	}
	
	inline void AbstractData::setBlock(const Block & b)
	{
		block=b;
	}
		
}
#endif

