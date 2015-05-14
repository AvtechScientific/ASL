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


#ifndef ASLFDADVECTIONDIFFUSION_H
#define ASLFDADVECTIONDIFFUSION_H

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
	
	/// Numerical method which computes multicomponent transport processes
	/**
		 \ingroup TransportProcesses
		 \ingroup NumMethods
		 
		 \f[ \partial_t c_i= D_i \Delta c_i - \nabla (\vec v c_i)
		    					-\nabla\left(\frac{c_i q}{k} \nabla(\phi+f_2)  \right)\f]
		 where
		 \param cData corresponds to \f$c_i\f$
		 \param diffusionCoefficient corresponds to \f$D_i\f$
		 \param efFactor1 corresponds to \f$k\f$
		 \param efFactor2 corresponds to \f$f_2\f$
		 \param efPhi corresponds to \f$\phi\f$
		 \param efChargeAnd corresponds to \f$q\f$	 
		 \param velocity corresponds to \f$\vec v\f$	 
	*/
	class FDAdvectionDiffusion: public SingleKernelNM
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPAbstractDataWithGhostNodes Field;
			
		private:	
			std::vector<Data> cData;
			std::vector<Data> cInternalData;			

			bool electricField;
			Field efPhi;
			Field efFactor1;
			Field efFactor2;
			std::vector<Field> efChargeAnd;

			Field velocity;
			/// advection process described by distribution function instead of velocity field
			Field distributionFunction;
			bool compressibilityCorrectionFlag;

			const VectorTemplate* vectorTemplate;

			int t;
			std::vector<acl::VectorOfElements> diffusionCoefficient;

			virtual void init0();
			virtual void postProcessing();
		public:			
			FDAdvectionDiffusion();
			FDAdvectionDiffusion(Data c, 
			                     const acl::VectorOfElements & dC, 
			                     const VectorTemplate* vT);
			void setDiffusionCoefficient(acl::VectorOfElements d, 
			                             unsigned int i = 0);
			inline const acl::VectorOfElements & getDiffusionCoefficient(unsigned int i=0) const;
			void setVectorTemplate(VectorTemplate* vT);
			inline const VectorTemplate* getVectorTemplate() const;
			void setElectricFieldParameters(Field phi,
			                                Field f1,
			                                Field f2,
			                                Field qAnd);
			void setVelocity(Field v, bool compressibilityCorrection=false);
			void setDistributionFunction(Field f);
			
			inline Field getVelocity();
			inline Field getDistributionFunction();
			inline std::vector<Data> & getData();
			void addComponent(Data c, acl::VectorOfElements & dC);
			void addComponent(Data c, acl::VectorOfElements & dC, Field qAnd);			
			
	};

	typedef std::shared_ptr<FDAdvectionDiffusion> SPFDAdvectionDiffusion;

	/**
		 \ingroup TransportProcesses
		 \ingroup NumMethods
		 
		 \f[ \partial_t c_i= D_i \Delta c_i - \nabla (\vec v c_i)\f]
		 where
		 \param cData corresponds to \f$c_i\f$
		 \param diffusionCoefficient corresponds to \f$D_i\f$
		 \param velocity corresponds to \f$\vec v\f$	 
	*/
	SPFDAdvectionDiffusion generateFDAdvectionDiffusion(SPDataWithGhostNodesACLData c, 
	                                                  double diffustionCoeff,
	                                                  SPAbstractDataWithGhostNodes v, 
	                                                  const VectorTemplate* vt,
	                                                  bool compressibilityCorrection = false);
	
	/**
		 \ingroup TransportProcesses
		 \ingroup NumMethods
		 
		 \f[ \partial_t c_i= D_i \Delta c_i \f]
		 where
		 \param cData corresponds to \f$ c_i \f$
		 \param diffusionCoefficient corresponds to \f$ D_i \f$
	*/
	SPFDAdvectionDiffusion generateFDAdvectionDiffusion(SPDataWithGhostNodesACLData c, 
	                                                    double diffustionCoeff,
	                                                    const VectorTemplate* vt);
	
// ------------------------- Implementation ------------------------

	inline FDAdvectionDiffusion::Field FDAdvectionDiffusion::getVelocity()
	{
		return velocity;
	}

	inline FDAdvectionDiffusion::Field FDAdvectionDiffusion::getDistributionFunction()
	{
		return distributionFunction;
	}
	
	inline std::vector<FDAdvectionDiffusion::Data> & FDAdvectionDiffusion::getData()
	{
		return cData;
	}
	
	inline const VectorTemplate* FDAdvectionDiffusion::getVectorTemplate() const
	{
		return vectorTemplate;
	}

	inline const acl::VectorOfElements & 
		FDAdvectionDiffusion::getDiffusionCoefficient(unsigned int i) const
	{
		return diffusionCoefficient[i];
	}
		
} // asl
#endif // ASLFDADVECTIONDIFFUSION_H
