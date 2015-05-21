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


#include "aslFDAdvectionDiffusionInhomogeneous.h"
#include "aslGenerators.h"
#include <acl/aclGenerators.h>
#include "acl/acl.h"
#include "math/aslTemplateVE.h"

#include <data/aslDataWithGhostNodes.h>
#include <acl/Kernels/aclKernel.h>
#include <math/aslTemplates.h>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>

namespace asl
{
	FDAdvectionDiffusionInhomogeneous::FDAdvectionDiffusionInhomogeneous():
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		cData(0),
		cInternalData(0),
		vectorTemplate(NULL),		
		diffusivity(0)
	{
	}

	FDAdvectionDiffusionInhomogeneous::FDAdvectionDiffusionInhomogeneous(Data c, ScalarField dC, const VectorTemplate* vT):
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		cData(1u,c),
		cInternalData(0u),
		vectorTemplate(vT),
		diffusivity(1u,dC)
	{
	}
	
	void FDAdvectionDiffusionInhomogeneous::init(){
		if (cData.size() != diffusivity.size())
			errorMessage("FDAdvectionDiffusion::init: some of compenents are underdefined");

		acl::TypeID type(getElementType(cData[0]->getDContainer()));
		
		cInternalData.resize(cData.size());
		for (unsigned int i(0); i < cData.size(); ++i)
			cInternalData[i] =clone(cData[i]);
		acl::VectorOfElements cnew(acl::generateVEPrivateVariable(1,type)); 
		for (unsigned int i(0); i < cData.size(); ++i){
			TemplateVE cT(*cData[i],d2q5());
			TemplateVE diffusivityT(*diffusivity[i],d2q5());
			(*kernel)<<cT.initValues;
			(*kernel)<<diffusivityT.initValues;
			(*kernel)<<(cnew=
						divAgradB(diffusivityT,cT)+
						cT.getValue(0));
			(*kernel)<<(assignmentSafe (cInternalData[i]->getSubContainer(),cnew));
		}		
		kernel->setup();
	}


	void FDAdvectionDiffusionInhomogeneous::execute()
	{
		kernel->compute();
		for (unsigned int i(0); i < cData.size(); ++i)
			swapBuffers(cData[i]->getDContainer(),cInternalData[i]->getDContainer());
	}
	
	void FDAdvectionDiffusionInhomogeneous::addComponent (Data c, ScalarField dC)
	{
		diffusivity.push_back(dC);
		cData.push_back (c);
	}
} //asl
