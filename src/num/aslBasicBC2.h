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


#ifndef ASLBASICBC2_H
#define ASLBASICBC2_H

#include "aslBCond.h"
#include <data/aslDataWithGhostNodes.h>
#include <acl/aclMath/aclVectorOfElementsDef.h>


namespace acl{
	class Kernel;
	typedef std::shared_ptr<Kernel> SPKernel;
}



namespace asl
{
	template <typename T> class UValue;
	class PositionFunction;
	typedef std::shared_ptr<PositionFunction> SPPositionFunction;



	/// Bondary condition that puts fixed value in each point
	/**
		\ingroup GenericBC		 
	*/
/*	class BCConstantValueMap:public BCondWithMap
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernel;			
		protected:
			Data data;
			acl::VectorOfElements value; 
		public:
			BCConstantValueMap(Data d, 
			                   const acl::VectorOfElements & v, 
			                   Data map);
			~BCConstantValueMap();
			virtual void execute();
			virtual void init();		
			void setValue(const acl::VectorOfElements & v);
	};
*/

	/// Bondary condition that puts fixed value in each point
	/**
		\ingroup GenericBC		 
	*/
/*	class BCValuePFMap:public BCondWithMap
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernel;			
		protected:
			Data data;
			SPPositionFunction value;
		public:
			BCValuePFMap(Data d, 
		  	             SPPositionFunction val, 
			             Data map);
			~BCValuePFMap();
			virtual void execute();
			virtual void init();		
			void setValue(SPPositionFunction v);
	};
*/	
	
	/// Bondary condition that makes fixed gradient, second order accuracy \ingroup GenericBC		 
	class BCConstantGradientMap2:public BCondWithMap
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernelCN;
			acl::SPKernel kernelGN;
		protected:
			Data  data;
			acl::VectorOfElements value; 
		public:
			BCConstantGradientMap2(Data d, 
			                       const acl::VectorOfElements & v, 
			                       Data map, 
			                       const VectorTemplate *const t);
			BCConstantGradientMap2(Data d, 
			                       const acl::VectorOfElements & v, 
			                       Data map, 
			                       Data computationalDomain, 
			                       const VectorTemplate *const t);
			~BCConstantGradientMap2();
			virtual void execute();
			virtual void init();			
			void setValue(const acl::VectorOfElements & v);	
	};
		
} //asl

#endif //ASLBASICBC2_H
