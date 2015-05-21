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


#ifndef ASLDATAWITHGHOSTNODES_H
#define ASLDATAWITHGHOSTNODES_H

#include "aslDataWrapper.h"
#include <iostream>
#include <fstream>

using  namespace std;

namespace asl
{

	///		 \ingroup DataFields
	///	\todo resolve consistency issue with AbstractData::setBlock()
	class AbstractDataWithGhostNodes: public AbstractData
	{
		protected:
			Block internalBlock;
			/// num of cells within the border
			unsigned int ghostBorder;
			/// flag represents whther the data acces on the borders is allowed
			bool bordersDataAcces;
			inline AbstractDataWithGhostNodes();					
			inline AbstractDataWithGhostNodes(const Block & b, int nGN = 1, bool bDA = true);			
		public:
			inline int getSubContainerOffset() const;	
			inline unsigned int getSubContainerSize() const;	
			inline const Block & getInternalBlock() const;
			inline const unsigned int getGhostBorder() const;
			virtual acl::VectorOfElements getSubContainer() = 0;	
	};

	///	 \ingroup DataFields
	template <typename V> class DataWithGhostNodes: public AbstractDataWithGhostNodes
	{
		private: 
			DataWrapper<V> dw;
		public:
			inline DataWithGhostNodes();
			inline DataWithGhostNodes(const Block & b, int nGN = 1, bool bDA = true);
			inline DataWithGhostNodes(DataWrapper<V> & d, int nGN = 1);
			virtual acl::VectorOfElements getSubContainer();	
			virtual const acl::VectorOfElements getEContainer();	
			virtual const acl::VectorOfElementsData getDContainer() const;				
			inline void setContainer(const V & cont);
			inline V & getContainer();	
	};

	typedef DataWithGhostNodes<acl::VectorOfElements> DataWithGhostNodesACL;
	typedef DataWithGhostNodes<acl::VectorOfElementsData> DataWithGhostNodesACLData;
//	typedef DataWithGhostNodes<vector<MemElement> > DataWithGhostNodesMem;

	typedef shared_ptr<DataWithGhostNodesACL> SPDataWithGhostNodesACL;
	typedef shared_ptr<DataWithGhostNodesACLData> SPDataWithGhostNodesACLData;
	typedef shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;

	///Creates new DataWithGhostNodesACLData with same data structure like \p d
	/// \ingroup DataFields
	SPDataWithGhostNodesACLData clone(SPDataWithGhostNodesACLData d);

	///Creates new DataWithGhostNodesACLData with same data structure like \p d containing \p n first elements
	/// \ingroup DataFields
	SPDataWithGhostNodesACLData clone(SPDataWithGhostNodesACLData d, unsigned int n);
	
	template <typename V> inline std::shared_ptr<DataWithGhostNodes<V>> 
		resizeGhostNodes(std::shared_ptr<DataWithGhostNodes<V>> d, unsigned int newGN);

//--------------------------- Implementation --------------------------

	AbstractDataWithGhostNodes::AbstractDataWithGhostNodes()
	{
	}

	AbstractDataWithGhostNodes::AbstractDataWithGhostNodes(const Block & b, int nGN,bool bDA): 
		AbstractData (offset(b, nGN)),
		internalBlock(b),
		ghostBorder(nGN),
		bordersDataAcces(bDA)
	{		
	}
		

	template <typename V> DataWithGhostNodes<V>::DataWithGhostNodes()
	{
	}
	
	template <typename V> DataWithGhostNodes<V>::DataWithGhostNodes(const Block & b, int nGN,bool bDA): 
		AbstractDataWithGhostNodes (b,nGN,bDA),
		dw(offset(b, nGN))
	{		
	}

	template <typename V> DataWithGhostNodes<V>::DataWithGhostNodes(DataWrapper<V> & d, int nGN): 
		AbstractDataWithGhostNodes (offset(d.getBlock(), -nGN), nGN, true),
		dw(d)
	{		
	}
		
	const Block & AbstractDataWithGhostNodes::getInternalBlock() const
	{
		return internalBlock;
	}

	const unsigned int AbstractDataWithGhostNodes::getGhostBorder() const
	{
		return ghostBorder;
	}

	inline int AbstractDataWithGhostNodes::getSubContainerOffset() const
	{
		Block b(offset(internalBlock, ghostBorder));		
		int nD(b.getSize().getSize());
		return b.c2i(AVec<int>(nD,ghostBorder));
	}
	
	inline unsigned int AbstractDataWithGhostNodes::getSubContainerSize() const
	{
		const Block b(offset(internalBlock, ghostBorder));
		int nD(b.getSize().getSize());
		int s(b.c2i(b.getSize()-AVec<int>(nD,ghostBorder+1)) - b.c2i(AVec<int>(nD,ghostBorder))+1);
		return s>0?s:0;
	}

	template <class V > inline void DataWithGhostNodes<V>::setContainer(const V & cont)
	{
		dw.setContainer (cont);
	}

	template <class V > inline V & DataWithGhostNodes<V>::getContainer()
	{
		return dw.getContainer();
	}
		
	template <typename V> inline std::shared_ptr<DataWithGhostNodes<V>> 
		resizeGhostNodes(std::shared_ptr<DataWithGhostNodes<V>> d, unsigned int newGN)
	{
		Block b(offset(d->getBlock(),newGN));
		std::shared_ptr<DataWithGhostNodes<V>> 
			nd(new DataWithGhostNodes<V>(b, 0));			
		nd->setContainer(d->getContainer());
		return nd;
	}

		
} //asl

#endif // ASLDATAWITHGHOSTNODES_H

