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


#include "aslFDAdvectionDiffusion.h"
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
	FDAdvectionDiffusion::FDAdvectionDiffusion():
		SingleKernelNM(acl::KERNEL_SIMDUA),
		cData(0),
		cInternalData(0),
		electricField(false),
		efChargeAnd(0),
		compressibilityCorrectionFlag(false),
		vectorTemplate(NULL)
	{
	}


	FDAdvectionDiffusion::FDAdvectionDiffusion(Data c, 
	                                           const acl::VectorOfElements & dC, 
	                                           const VectorTemplate* vT):
		SingleKernelNM(acl::KERNEL_SIMDUA),
		cData(1u, c),
		cInternalData(0u),
		electricField(false),		
		efChargeAnd(0u),
		compressibilityCorrectionFlag(false),
		vectorTemplate(vT),
		diffusionCoefficient(1u, dC)
	{
	}

		
	void FDAdvectionDiffusion::setElectricFieldParameters(Field phi, Field f1, Field f2, Field qAnd)
	{
		electricField = true;
		efPhi = phi;
		efFactor1 = f1;
		efFactor2 = f2;
		efChargeAnd.resize(1);
		efChargeAnd[0] = qAnd;
	}

		
	void FDAdvectionDiffusion::init0()
	{
		if (electricField)
		{
			if (cData.size() != diffusionCoefficient.size() || cData.size() != efChargeAnd.size())
				errorMessage("FDAdvectionDiffusion::init - some of compenents are underdefined");
		}
		else
		{
			if (cData.size() != diffusionCoefficient.size())
				errorMessage("FDAdvectionDiffusion::init - some of compenents are underdefined");
		}

		acl::TypeID type(getElementType(cData[0]->getDContainer()));
		
		cInternalData.resize(cData.size());
		for (unsigned int i(0); i < cData.size(); ++i)
			cInternalData[i] = clone(cData[i]);
		acl::VectorOfElements cnew(acl::generateVEPrivateVariable(1, type)); 
		std::unique_ptr<TemplateVE> f1;
		std::unique_ptr<TemplateVE> f2;
		std::unique_ptr<TemplateVE> phi;


		if (electricField)
		{
			f1.reset(new TemplateVE(*efFactor1, *vectorTemplate));
			f2.reset(new TemplateVE(*efFactor2, *vectorTemplate));				
			phi.reset(new TemplateVE(*efPhi, *vectorTemplate));				
			(*kernel) << f1->initValues;
			(*kernel) << f2->initValues;
			(*kernel) << phi->initValues;
		}

		unsigned int nd(nD(*vectorTemplate));

		vector<TemplateVE> velocityT(nd);
		if(velocity.get() != 0)
			for(unsigned int i(0); i < nd; ++i)
			{
				velocityT[i].init(*velocity, *vectorTemplate, i);
				*(kernel) << velocityT[i].initValues;
			}
		
		TemplateVE cT;
		for (unsigned int i(0); i < cData.size(); ++i)
		{
			cT.init(*cData[i], *vectorTemplate);
			(*kernel) << cT.initValues;
			(*kernel) << (cnew = 
			              diffusionCoefficient[i] * laplas(cT) + 
			              cT.getValue(0));
			if(velocity.get() != 0)
			{
				(*kernel) << (cnew -= divProduct(velocityT, cT));
				if(compressibilityCorrectionFlag)
				{
					(*kernel) << (cnew /= (1. - div(velocityT)));
					// improovment of stability
					if(nd == 3)
					{
						auto vl(acl::cat(subVE(velocityT[0].values,0), 
						                 subVE(velocityT[1].values,0), 
						                 subVE(velocityT[2].values,0)));
						(*kernel) << (cnew += laplas(cT) * (vl*vl)*.5);
					}
				}
			}
//			if(distributionFunction.get() != 0)
//			{
//			}
			if (electricField)
			{
				TemplateVE q(*efChargeAnd[i], *vectorTemplate);
				(*kernel) << q.initValues;
				(*kernel) << (cnew -= divAgradB(cT*q/(*f1), (*phi)/*+(*f2)*/)); 
			}
			(*kernel) << (assignmentSafe(cInternalData[i]->getSubContainer(), cnew));
		}		
	}

		
	void FDAdvectionDiffusion::postProcessing()
	{
		for (unsigned int i(0); i < cData.size(); ++i)
			swapBuffers(cData[i]->getDContainer(), cInternalData[i]->getDContainer());
	}
		

	void FDAdvectionDiffusion::setDiffusionCoefficient(acl::VectorOfElements dC, 
	                                                   unsigned int i)
	{
		copy(dC,diffusionCoefficient[i]);
	}

		
	void FDAdvectionDiffusion::addComponent(Data c, acl::VectorOfElements & dC)
	{
		if (electricField)
			errorMessage("FDAdvectionDiffusion::addComponent: \"Electric field\" is swiched on \n therefore \"qAnd\" value should be specified");
		unsigned int n(diffusionCoefficient.size());
		diffusionCoefficient.resize(n+1);
		setDiffusionCoefficient(dC,n);
		
		cData.push_back(c);
	}

		
	void FDAdvectionDiffusion::addComponent(Data c, 
	                                        acl::VectorOfElements & dC,
	                                        Field qAnd)
	{
		unsigned int n(diffusionCoefficient.size());
		diffusionCoefficient.resize(n+1);
		setDiffusionCoefficient(dC,n);
		cData.push_back(c);
		efChargeAnd.push_back(qAnd);
	}

	void FDAdvectionDiffusion::setVelocity(Field v, bool cCF)
	{
		velocity=v;
		compressibilityCorrectionFlag=cCF;
	}

	void FDAdvectionDiffusion::setDistributionFunction(Field f)
	{
		if(f->getEContainer().size() != vectorTemplate->vectors.size())
			errorMessage("FDAdvectionDiffusion::setDistributionFunction: the distrubution function has wrong number of components");
		distributionFunction=f;
	}

		SPFDAdvectionDiffusion generateFDAdvectionDiffusion(SPDataWithGhostNodesACLData c,
		                                                    double diffusionCoeff,            
		                                                    SPAbstractDataWithGhostNodes v, 
		                                                    const VectorTemplate* vt, 
		                                                    bool compressibilityCorrection)
	{
		auto nm(make_shared<FDAdvectionDiffusion> (c, acl::generateVEConstant(diffusionCoeff), vt));
		nm->setVelocity(v,compressibilityCorrection);
		return nm;
	}

		SPFDAdvectionDiffusion generateFDAdvectionDiffusion(SPDataWithGhostNodesACLData c,
		                                                    double diffusionCoeff, 
		                                                    const VectorTemplate* vt)
	{
		auto nm(make_shared<FDAdvectionDiffusion> (c, acl::generateVEConstant(diffusionCoeff), vt));
		return nm;
	}
		
		
} //asl
