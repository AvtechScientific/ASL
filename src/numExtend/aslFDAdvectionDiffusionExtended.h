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

#include "../num/aslNumMethod.h"


namespace acl{
	class Kernel;
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
	
	/// Numerical method which computes multicomponent transport processes
	/**
		 \ingroup TransportProcesses
		 \f[ \partial_t c_i= D_i \Delta c_i 
		    					-\nabla\left(\frac{c_i q}{k} \nabla(\phi+f_2)  \right)\f]
		 where
		 \param cData corresponds to \f$c_i\f$
		 \param diffusionCoefficient corresponds to \f$D_i\f$
		 \param efFactor1 corresponds to \f$k\f$
		 \param efFactor2 corresponds to \f$f_2\f$
		 \param efPhi corresponds to \f$\phi\f$
		 \param efChargeAnd corresponds to \f$q\f$

		 This class contains filtering which forbits the local value of 
		 \f$ \frac{q}{k} \nabla(\phi+f_2) \f$ to be larger than \f$0.2\f$
	 */
	class FDAdvectionDiffusionExtended: public NumMethod
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPAbstractDataWithGhostNodes ScalarField;
		private:
			std::unique_ptr<acl::Kernel> kernel; 
			
			std::vector<Data> cData;
			std::vector<Data> cInternalData;			

			bool electricField;
			ScalarField efPhi;
			ScalarField efFactor1;
			ScalarField efFactor2;
			std::vector<ScalarField> efChargeAnd;
				
			const VectorTemplate* vectorTemplate;

			int t;
			std::vector<double> diffusionCoefficient;
		public:			
			FDAdvectionDiffusionExtended();
			FDAdvectionDiffusionExtended(Data c, double dC,const VectorTemplate* vT);
			void setDiffusionCoefficient(double d, unsigned int i=0);
			double getDiffusionCoefficient(unsigned int i=0);
			void setVectorTemplate(VectorTemplate* vT);
			VectorTemplate setVectorTemplate(const VectorTemplate* vT);
			void setElectricFieldParameters(ScalarField phi, ScalarField f1, ScalarField f2, ScalarField qAnd);
			virtual void init();
			virtual void execute();
			void addComponent(Data c, double dC);
			void addComponent(Data c, double dC, ScalarField qAnd);			
			
	};

	typedef std::shared_ptr<FDAdvectionDiffusionExtended> SPFDAdvectionDiffusionExtended;

		
} // asl
#endif // ASLFDADVECTIONDIFFUSION_H
