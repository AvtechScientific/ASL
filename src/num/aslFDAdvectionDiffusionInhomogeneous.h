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


#ifndef ASLFDADVECTIONDIFFUSIONINHOMOGENEOUS_H
#define ASLFDADVECTIONDIFFUSIONINHOMOGENEOUS_H

#include "aslNumMethod.h"


namespace acl
{
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
		 \f[ \partial_t [cData]_i= \nabla_\alpha [diffusivity]_i \nabla_\alpha [cData]_i  \f]
	 */
	class FDAdvectionDiffusionInhomogeneous: public NumMethod
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPAbstractDataWithGhostNodes ScalarField;
		private:
			std::unique_ptr<acl::Kernel> kernel; 
			
			std::vector<Data> cData;
			std::vector<Data> cInternalData;			
				
			const VectorTemplate* vectorTemplate;

			int t;
			std::vector<ScalarField> diffusivity;
		public:			
			FDAdvectionDiffusionInhomogeneous();
			FDAdvectionDiffusionInhomogeneous(Data c, ScalarField dC, const VectorTemplate* vT);
			void setVectorTemplate(VectorTemplate* vT);
			VectorTemplate setVectorTemplate(const VectorTemplate* vT);
			virtual void init();
			virtual void execute();
			void addComponent(Data c, ScalarField dC);
	};

	typedef std::shared_ptr<FDAdvectionDiffusionInhomogeneous> SPFDAdvectionDiffusionInhomogeneous;

		
} // asl
#endif // ASLFDADVECTIONDIFFUSIONINHOMOGENEOUS_H
