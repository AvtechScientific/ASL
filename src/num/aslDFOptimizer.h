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


#ifndef ASLDFOPTIMIZER_H
#define ASLDFOPTIMIZER_H

#include "aslSingleKernelNM.h"

namespace acl
{
	class VectorOfElementsData;
	class VectorOfElements;
}

namespace asl
{
	class VectorTemplate;
	template <typename V> class DataWithGhostNodes;
	typedef DataWithGhostNodes<acl::VectorOfElementsData> DataWithGhostNodesACLData;
	typedef std::shared_ptr<DataWithGhostNodesACLData> SPDataWithGhostNodesACLData;
	class AbstractDataWithGhostNodes;
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;
	
	/// Numerical method which makes changes in the input map and produces map suitbale and optimal for use in BC
	/**
		 \ingroup LevelSet
		 \ingroup NumMethods
	 
	*/
	class DFOptimizer: public SingleKernelNM
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPAbstractDataWithGhostNodes Field;
		private:	
			Data inData;

			const VectorTemplate* vectorTemplate;

			virtual void init0();
		public:			
			DFOptimizer();
			DFOptimizer(Data inD, const VectorTemplate* vT);
			void setVectorTemplate(VectorTemplate* vT);
			inline const VectorTemplate* getVectorTemplate() const;
			inline Data & getInData();
	};

	typedef std::shared_ptr<DFOptimizer> SPDFOptimizer;

	void optimizeMap(SPDataWithGhostNodesACLData c, const VectorTemplate* vt);
	
// ------------------------- Implementation ------------------------
	
	inline DFOptimizer::Data & DFOptimizer::getInData()
	{
		return inData;
	}
	
	inline const VectorTemplate* DFOptimizer::getVectorTemplate() const
	{
		return vectorTemplate;
	}
		
} // asl
#endif // ASLFDADVECTIONDIFFUSION_H
