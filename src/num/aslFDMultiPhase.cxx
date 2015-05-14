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


#include "aslFDMultiPhase.h"
#include <aslGenerators.h>
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include <math/aslTemplateVE.h>
#include <acl/aclMath/aclVectorOfElements.h>

#include <data/aslDataWithGhostNodes.h>
#include <acl/Kernels/aclKernel.h>
#include <math/aslTemplates.h>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>

namespace asl
{
	FDMultiPhase::FDMultiPhase():
		SingleKernelNM(acl::KERNEL_SIMDUA),
		cData(0),
		cInternalData(0),
		compressibilityCorrectionFlag(false),
		vectorTemplate(NULL)
	{
	}


	FDMultiPhase::FDMultiPhase(Data c, 
	                           const acl::VectorOfElements & dC,
	                           const acl::VectorOfElements & rC,
	                           const VectorTemplate* vT):
		SingleKernelNM(acl::KERNEL_SIMDUA),
		cData(1u, c),
		cInternalData(0u),
		compressibilityCorrectionFlag(false),
		vectorTemplate(vT),
		diffusionCoefficient(dC),
		repulsionConstant(rC)
	{
	}
		
	void FDMultiPhase::init0()
	{

		acl::TypeID type(getElementType(cData[0]->getDContainer()));
		
		cInternalData.resize(cData.size());
		for (unsigned int i(0); i < cData.size(); ++i)
			cInternalData[i] = clone(cData[i]);
		acl::VectorOfElements cnew(acl::generateVEPrivateVariable(1, type)); 
		std::unique_ptr<TemplateVE> f1;
		std::unique_ptr<TemplateVE> f2;
		std::unique_ptr<TemplateVE> phi;

		unsigned int nd(nD(*vectorTemplate));

		vector<TemplateVE> velocityT(nd);
		if(velocity.get() != 0)
			for(unsigned int i(0); i < nd; ++i)
			{
				velocityT[i].init(*velocity, *vectorTemplate, i);
				*(kernel) << velocityT[i].initValues;
			}
		
		TemplateVE cT;
		auto vOne(acl::generateVEConstantN(vectorTemplate->vectors.size(),1.));
		for (unsigned int i(0); i < cData.size(); ++i)
		{
			cT.init(*cData[i], *vectorTemplate);
			(*kernel) << cT.initValues;
			TemplateVE cTAbove(productOfElements((cT.values-vOne),.5*(sign(cT.values-vOne)+vOne)),*vectorTemplate);
			TemplateVE cTBelow(productOfElements(cT.values,.5*(vOne-sign(cT.values))),*vectorTemplate);			
			
			(*kernel) << (cnew = 
			              diffusionCoefficient * laplas(cT) +
			              0.15 * laplas(cTAbove) +
			              0.15 * laplas(cTBelow) -                            
			              repulsionConstant* laplas(cT - cT * cT) +			              
			              cT.getValue(0));
			if(velocity.get() != 0)
			{
				(*kernel) << (cnew -= divProduct(velocityT, cT));
				if(compressibilityCorrectionFlag)
				{
					(*kernel) << (cnew /= (1. - div(velocityT)));
				}
			}
			(*kernel) << (assignmentSafe(cInternalData[i]->getSubContainer(), cnew));
		}		
	}

		
	void FDMultiPhase::postProcessing()
	{
		for (unsigned int i(0); i < cData.size(); ++i)
			swapBuffers(cData[i]->getDContainer(), cInternalData[i]->getDContainer());
	}
		

	void FDMultiPhase::setDiffusionCoefficient(acl::VectorOfElements dC)
	{
		copy(dC,diffusionCoefficient);
	}

		
	void FDMultiPhase::addComponent(Data c)
	{	
		cData.push_back(c);
	}

		
	void FDMultiPhase::setVelocity(Field v, bool cCF)
	{
		velocity=v;
		compressibilityCorrectionFlag=cCF;
	}

	SPFDMultiPhase generateFDMultiPhase(SPDataWithGhostNodesACLData c,
	                                    SPAbstractDataWithGhostNodes v, 
	                                    const VectorTemplate* vt, 
	                                    bool compressibilityCorrection)
	{
		auto nm(make_shared<FDMultiPhase> (c, 
		                                   acl::generateVEConstant(0.07),
		                                   acl::generateVEConstant(0.1),
		                                   vt));
		nm->setVelocity(v,compressibilityCorrection);
		return nm;
	}

	
		
} //asl
