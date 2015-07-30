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


#ifndef ASLDATARESAMPLING_H
#define ASLDATARESAMPLING_H

#include "aslSingleKernelNM.h"
#include "math/aslVectors.h"
//#include <CL/cl.hpp>
// Supply "cl.hpp" with ASL, since it is not present in OpenCL 2.0
// Remove the file after switching to OpenCL 2.1
#include "acl/cl.hpp"

namespace acl
{
	class VectorOfElementsData;
}

namespace asl
{
	class VectorTemplate;
	template <typename V> class DataWithGhostNodes;
	typedef DataWithGhostNodes<acl::VectorOfElementsData> DataWithGhostNodesACLData;
	typedef std::shared_ptr<DataWithGhostNodesACLData> SPDataWithGhostNodesACLData;
	class AbstractDataWithGhostNodes;
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;
	
	/// Algorithm for generation of coarsed dataset
	/**
		 \ingroup NumMethods		 
		 
		 \todo make and test
	*/
	class DataCoarser: public SingleKernelNM
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			
		private:	
			Data dataIn;
			Data dataOut;
			const VectorTemplate* vectorTemplate;
	
			virtual void init0();
		public:			
			DataCoarser();
			DataCoarser(Data dIn);
			inline Data getDataOut();
	};

	typedef std::shared_ptr<DataCoarser> SPDataCoarser;

	/// \ingroup DataUtilities
	SPDataWithGhostNodesACLData coarseData(SPDataWithGhostNodesACLData d);

	/// Algorithm for generation of coarsed dataset
	/**
		 \ingroup NumMethods		 
		 
		 \todo make and test
	*/
	class DataClipper: public SingleKernelNM
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			
		private:	
			Data dataIn;
			Data dataOut;
			AVec<int> a0;
			AVec<int> aE;
	
			virtual void init0();
		public:			
			DataClipper();
			DataClipper(Data dIn, AVec<int> a0, AVec<int> aE);
			inline Data getDataOut();
	};
	
	typedef std::shared_ptr<DataClipper> SPDataClipper;

	/// \ingroup DataUtilities
	inline SPDataWithGhostNodesACLData clipData(SPDataWithGhostNodesACLData d,
	                                             AVec<int> a0, 
	                                             AVec<int> aE);

//----------------------------- Implementation -----------------------

	DataCoarser::Data DataCoarser::getDataOut()
	{
		return dataOut;
	}

	DataClipper::Data DataClipper::getDataOut()
	{
		return dataOut;
	}

	inline SPDataWithGhostNodesACLData coarseData(SPDataWithGhostNodesACLData d)
	{
		DataCoarser dc(d);
		dc.init();
		dc.execute();
		return dc.getDataOut();
		
	}

	inline SPDataWithGhostNodesACLData clipData(SPDataWithGhostNodesACLData d,
	                                     AVec<int> a0, 
	                                     AVec<int> aE)
	{
		asl::DataClipper dcl(d, a0,aE);
		dcl.init();
		dcl.execute();
		return dcl.getDataOut();
	}

	
} // asl
#endif // ASLDATARESAMPLING_H
