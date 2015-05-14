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


#include "aslFDAdvectionDiffusion2.h"
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
	FDAdvectionDiffusion2::FDAdvectionDiffusion2():
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		cData(0),
		cInternalData(0),
		electricField(false),
		efChargeAnd(0),
		vectorTemplate(NULL),		
		diffusionCoefficient(0)
	{
	}

	FDAdvectionDiffusion2::FDAdvectionDiffusion2(Data c, double dC, const VectorTemplate* vT):
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		cData(1u,c),
		cInternalData(0u),
		electricField(false),		
		efChargeAnd(0u),
		vectorTemplate(vT),
		diffusionCoefficient(1u,dC)
	{
	}


	void FDAdvectionDiffusion2::setElectricFieldParameters(ScalarField phi, ScalarField f1, ScalarField f2, ScalarField qAnd)
	{
		electricField=true;
		efPhi=phi;
		efFactor1=f1;
		efFactor2=f2;
		efChargeAnd.resize(1);
		efChargeAnd[0] =qAnd;
	}
	
	void FDAdvectionDiffusion2::init(){
		if (electricField)
		{
			if (cData.size() != diffusionCoefficient.size() || cData.size() != efChargeAnd.size())
				errorMessage("FDAdvectionDiffusion2::init: some of compenents are underdefined");
		}
		else
		{
			if (cData.size() != diffusionCoefficient.size())
				errorMessage("FDAdvectionDiffusion2::init: some of compenents are underdefined");
		}

		acl::TypeID type(getElementType(cData[0]->getDContainer()));
		
		cInternalData.resize(cData.size());
		for (unsigned int i(0); i < cData.size(); ++i)
			cInternalData[i] =clone(cData[i]);
		acl::VectorOfElements cnew(acl::generateVEPrivateVariable(1,type)); 
		std::unique_ptr<TemplateVE> f1;
		std::unique_ptr<TemplateVE> f2;
		std::unique_ptr<TemplateVE> phi;
		if (electricField)
		{
			f1.reset(new TemplateVE(*efFactor1,d2q5()));
			f2.reset(new TemplateVE(*efFactor2,d2q5()));				
			phi.reset(new TemplateVE(*efPhi,d2q5()));				
			(*kernel)<<f1->initValues;
			(*kernel)<<f2->initValues;
			(*kernel)<<phi->initValues;
		}
		
		TemplateVE cT;
		for (unsigned int i(0); i < cData.size(); ++i){
			cT.init(*cData[i],d2q5());
			(*kernel)<<cT.initValues;
			(*kernel)<<(cnew=
						acl::generateVEVariableR(diffusionCoefficient[i])*laplas(cT));
			if (electricField)
			{
				TemplateVE q(*efChargeAnd[i],d2q5());
				(*kernel)<<q.initValues;
				(*kernel)<<(cnew-=divAgradB(cT*q/(*f1),(*phi)/*+(*f2)*/)); 
			}
			(*kernel)<<(assignmentSafe (cInternalData[i]->getSubContainer(),cInternalData[i]->getSubContainer()+acl::generateVEConstant(2.)*cnew));
		}		
		kernel->setup();
	}


	void FDAdvectionDiffusion2::execute()
	{
		kernel->compute();
		for (unsigned int i(0); i < cData.size(); ++i)
			swapBuffers(cData[i]->getDContainer(),cInternalData[i]->getDContainer());
	}


	void FDAdvectionDiffusion2::setDiffusionCoefficient (double dC, unsigned int i)
	{
		diffusionCoefficient[i] =dC;
	}
	
	void FDAdvectionDiffusion2::addComponent (Data c, double dC)
	{
		if (electricField)
			errorMessage("FDAdvectionDiffusion2::addComponent: \"Electric field\" is swiched on \n therefore \"qAnd\" value should be specified");
		diffusionCoefficient.push_back (dC);
		cData.push_back (c);
	}
	void FDAdvectionDiffusion2::addComponent (Data c, double dC,ScalarField qAnd)
	{
		diffusionCoefficient.push_back (dC);
		cData.push_back (c);
		efChargeAnd.push_back (qAnd);
	}
} //asl
