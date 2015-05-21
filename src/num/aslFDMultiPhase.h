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


#ifndef ASLFDMULTIPHASE_H
#define ASLFDMULTIPHASE_H

#include "aslSingleKernelNM.h"
#include "acl/aclMath/aclVectorOfElementsDef.h"

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
	
	/// Numerical method which computes multiphase transport processes
	/**
		 \ingroup TransportProcesses
		 \ingroup NumMethods
		 
		 \f[ \partial_t c_i= D \Delta c_i - \nabla (\vec v c_i)
		    					-\nabla\left( a c_i \sum_{j\neq i}\nabla  c_j  \right)\f]
		 where
		 \param cData corresponds to \f$c_i\f$
		 \param diffusionCoefficient corresponds to \f$D\f$
		 \param repulsionconstant is repulsion constant \f$a\f$ 
		 \param velocity corresponds to \f$\vec v\f$	 
	*/
	class FDMultiPhase: public SingleKernelNM
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPAbstractDataWithGhostNodes Field;
			
		private:	
			std::vector<Data> cData;
			std::vector<Data> cInternalData;			

			Field velocity;
			bool compressibilityCorrectionFlag;

			const VectorTemplate* vectorTemplate;

			int t;
			acl::VectorOfElements diffusionCoefficient;
			acl::VectorOfElements repulsionConstant;

			virtual void init0();
			virtual void postProcessing();
		public:			
			FDMultiPhase();
			FDMultiPhase(Data c, 
			             const acl::VectorOfElements & dC,
			             const acl::VectorOfElements & rC,
			             const VectorTemplate* vT);
			void setDiffusionCoefficient(acl::VectorOfElements d);
			inline const acl::VectorOfElements & getDiffusionCoefficient() const;
			void setRepulsionConstant(acl::VectorOfElements d);
			inline const acl::VectorOfElements & getRepulsionConstant() const;
			void setVectorTemplate(VectorTemplate* vT);
			inline const VectorTemplate* getVectorTemplate() const;
			void setVelocity(Field v, bool compressibilityCorrection=false);
			inline Field getVelocity();
			inline std::vector<Data> & getData();
			void addComponent(Data c);
	};

	typedef std::shared_ptr<FDMultiPhase> SPFDMultiPhase;

	/**
		 \ingroup TransportProcesses
		 \ingroup NumMethods
		 
		 \f[ \partial_t c_i= D \Delta c_i - \nabla (\vec v c_i) - 
			                 \nabla\left( a c_i \sum_{j\neq i}\nabla  c_j\f]
		 where
		 \param cData corresponds to \f$c_i\f$
		 \param diffusionCoefficient corresponds to \f$D_i\f$
		 \param velocity corresponds to \f$\vec v\f$	 
	*/
	SPFDMultiPhase generateFDMultiPhase(SPDataWithGhostNodesACLData c, 
	                                            SPAbstractDataWithGhostNodes v, 
	                                            const VectorTemplate* vt,
	                                            bool compressibilityCorrection = false);
	
	/**
		 \ingroup TransportProcesses
		 \ingroup NumMethods
		 
		 \f[ \partial_t c_i= D_i \Delta c_i \f]
		 where
		 \param cData corresponds to \f$ c_i \f$
			 
		 \f$ D_i is diffusionCoefficient \f$,  \f$a\f$ is repulsion constant			 
	*/
	SPFDMultiPhase generateFDAdvectionDiffusion(SPDataWithGhostNodesACLData c, 
	                                            double diffustionCoeff,
	                                            const VectorTemplate* vt);
	
// ------------------------- Implementation ------------------------

	inline FDMultiPhase::Field FDMultiPhase::getVelocity()
	{
		return velocity;
	}
	
	inline std::vector<FDMultiPhase::Data> & FDMultiPhase::getData()
	{
		return cData;
	}
	
	inline const VectorTemplate* FDMultiPhase::getVectorTemplate() const
	{
		return vectorTemplate;
	}

	inline const acl::VectorOfElements & 
		FDMultiPhase::getDiffusionCoefficient() const
	{
		return diffusionCoefficient;
	}
		
} // asl
#endif // ASLFDMULTIPHASE_H
