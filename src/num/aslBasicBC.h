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


#ifndef ASLBASICBC_H
#define ASLBASICBC_H

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
	class BCConstantValue:public BCond
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernel;			
		protected:
			Data data;
			acl::VectorOfElements value; 
		public:
			BCConstantValue(Data d, const acl::VectorOfElements & v);
			virtual void execute();
			virtual void init();		
			void setValue(const acl::VectorOfElements & v);
	};

	/// Bondary condition that puts fixed value in each point
	/**
		\ingroup GenericBC		 
	*/
	class BCConstantValueMap:public BCondWithMap
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

	/// Bondary condition that puts fixed value in each boundary point
	/**
		\ingroup GenericBC		 
	*/	
	class BCConstantValueMiddlePointMap:public BCondWithMap
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernel;			
		protected:
			Data data;
			acl::VectorOfElements value; 
		public:
			BCConstantValueMiddlePointMap(Data d, 
			                              const acl::VectorOfElements & v, 
			                              Data map,
			                              const VectorTemplate *const t);
			~BCConstantValueMiddlePointMap();
			virtual void execute();
			virtual void init();		
			void setValue(const acl::VectorOfElements & v);
	};

	/// Bondary condition that puts fixed value in each point
	/**
		\ingroup GenericBC		 
	*/
	class BCValuePFMap:public BCondWithMap
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
	
	
	/// Bondary condition that puts fixed value in each point
	/**
		\ingroup GenericBC		 
	*/
	SPBCond generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                double v, 
	                                const std::vector<SlicesNames> & sl);
	/// Bondary condition that puts fixed value in each point
	/**
		\ingroup GenericBC		 
	*/
	SPBCond generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                UValue<double> & v, 
	                                const std::vector<SlicesNames> & sl);
	/// Bondary condition that puts fixed value in each point
	/**
		\ingroup GenericBC		 
	*/
	SPBCond generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                UValue<AVec<float>> & v, 
	                                const std::vector<SlicesNames> & sl);
	/// Bondary condition that puts fixed value in each point
	/**
		\ingroup GenericBC		 
	*/
	SPBCond generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                AVec<> v, 
	                                const std::vector<SlicesNames> & sl);
	
	/// Bondary condition that puts fixed value in each point \ingroup GenericBC		 
	SPNumMethod generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                   double v, 
	                                   SPAbstractDataWithGhostNodes map);
	
	/// Bondary condition that puts fixed value in each point \ingroup GenericBC		 
	SPNumMethod generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                   AVec<> v, 
	                                   SPAbstractDataWithGhostNodes map);

	/// Bondary condition that puts fixed value in each point \ingroup GenericBC		 
	SPNumMethod generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                   SPPositionFunction v, 
	                                   SPAbstractDataWithGhostNodes map);

	/// Bondary condition that puts fixed value in each point \ingroup GenericBC		 
	SPNumMethod generateBCConstantValueMiddlePoint(SPAbstractDataWithGhostNodes d, 
	                                               double v, 
	                                               SPAbstractDataWithGhostNodes map,
	                                               const VectorTemplate *const t);

	/// Bondary condition that puts fixed value in each point \ingroup GenericBC		 
	SPNumMethod generateBCConstantValueMiddlePoint(SPAbstractDataWithGhostNodes d, 
	                                               AVec<> v, 
	                                               SPAbstractDataWithGhostNodes map,
	                                               const VectorTemplate *const t);
	

	/// Bondary condition that makes fixed gradient	\ingroup GenericBC		 
	class BCConstantGradient:public BCond
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernel;			
		protected:
			Data  data;
			acl::VectorOfElements value; 
		public:
			BCConstantGradient(Data d, 
			               const acl::VectorOfElements & v, 
			               const VectorTemplate *const t);
			virtual void execute();
			virtual void init();			
			void setValue(const acl::VectorOfElements & value);	
	};

	/// Bondary condition that makes fixed gradient	\ingroup GenericBC		 
	class BCConstantGradientMap:public BCondWithMap
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernel;			
		protected:
			Data  data;
			acl::VectorOfElements value; 
		public:
			BCConstantGradientMap(Data d, 
			                  const acl::VectorOfElements & v, 
			                  Data map, 
			                  const VectorTemplate *const t);
			BCConstantGradientMap(Data d, 
			                  const acl::VectorOfElements & v, 
			                  Data map, 
			                  Data computationalDomain, 
			                  const VectorTemplate *const t);
			~BCConstantGradientMap();
			virtual void execute();
			virtual void init();			
			void setValue(const acl::VectorOfElements & v);	
	};
	
	/// Bondary condition that makes fixed gradient	\ingroup GenericBC		 
	SPBCond generateBCConstantGradient(SPAbstractDataWithGhostNodes d, 
	                               double v,
	                               const VectorTemplate *const t,	                               
	                               const std::vector<SlicesNames> & sl);

	/// Bondary condition that makes fixed gradient	\ingroup GenericBC		 
	SPNumMethod generateBCConstantGradient(SPAbstractDataWithGhostNodes d, 
	                                       double v, 
	                                       SPAbstractDataWithGhostNodes map,
	                                       const VectorTemplate *const t);
	/// Bondary condition that makes fixed gradient	\ingroup GenericBC		 
	SPNumMethod generateBCConstantGradient(SPAbstractDataWithGhostNodes d, 
	                                       double v, 
	                                       SPAbstractDataWithGhostNodes map,
	                                       SPAbstractDataWithGhostNodes computatinalDomain,
	                                       const VectorTemplate *const t);

	/// Bondary condition that makes fixed gradient	\ingroup GenericBC		 
	SPNumMethod generateBCConstantGradient(SPAbstractDataWithGhostNodes d, 
	                                       AVec<> v, 
	                                       SPAbstractDataWithGhostNodes map,
	                                       const VectorTemplate *const t);

	/// Bondary condition that makes fixed gradient, second order accuracy \ingroup GenericBC		 
	SPNumMethod generateBCConstantGradient2(SPAbstractDataWithGhostNodes d, 
	                                        double v, 
	                                        SPAbstractDataWithGhostNodes map,
	                                        const VectorTemplate *const t);
	/// Bondary condition that makes fixed gradient, second order accuracy \ingroup GenericBC		 
	SPNumMethod generateBCConstantGradient2(SPAbstractDataWithGhostNodes d, 
	                                        double v, 
	                                        SPAbstractDataWithGhostNodes map,
	                                        SPAbstractDataWithGhostNodes computatinalDomain,
	                                        const VectorTemplate *const t);
	
	/// Bondary condition that makes fixed gradient, second order accuracy \ingroup GenericBC		 
	SPNumMethod generateBCConstantGradient2(SPAbstractDataWithGhostNodes d, 
	                                        AVec<> v, 
	                                        SPAbstractDataWithGhostNodes map,
	                                        const VectorTemplate *const t);
	
	/// Bondary condition that adds fixed value to one in each point
	/**
		\ingroup GenericBC	
	*/
	class BCConstantSource:public BCond
	{		
		public:
			typedef SPDataWithGhostNodesACLData Data;
			acl::SPKernel kernel;			
		protected:
			Data data;
			cl_double value; 
		public:
			BCConstantSource(Data d, cl_double v);
			virtual void execute();
			virtual void init();			
			void setValue(cl_double value);
	};

	
	/// Bondary condition that copies directly the values from one data to another
	/**
		\ingroup GenericBC		 
	*/
	class BCDirectCopier:public BCondConnector
	{		
		public:
			typedef SPDataWithGhostNodesACLData Data;
			acl::SPKernel kernel;			
		protected:
			Data & source;
			Data & destination;
		public:
			BCDirectCopier(Data dSource,Data dDestination);
			virtual void execute();
			virtual void init();			
	};

	/// Bondary condition that puts fixed value in each point uses Slices
	/**
		\ingroup GenericBC		 
	*/
	class BCSConstantValue:public BCondSlice
	{		
		public:
			typedef SPDataWithGhostNodesACLData Data;
		protected:
			acl::SPKernel kernel;			
			Data data;
			cl_double value; 
		public:
			BCSConstantValue(Data d, cl_double v);
			virtual void execute();
			virtual void init();			
	};

	
	
} //asl

#endif //ASLBASICBC_H
