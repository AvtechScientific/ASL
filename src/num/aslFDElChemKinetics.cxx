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


#include "aslFDElChemKinetics.h"
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

	FDBVKinetics::FDBVKinetics(Data a0,
	                           double n0,
	                           Data aI,
	                           double nI,
	                           Data phiS,
	                           Field phi,
	                           const Param & j0,
	                           const Param & b,
	                           double n):
		SingleKernelMapNM(acl::KERNEL_SIMD),
		kernelJ(make_shared<acl::Kernel>(acl::KERNEL_SIMD)),
//		SingleKernelMapNM(acl::KERNEL_BASIC),
//		kernelJ(make_shared<acl::Kernel>(acl::KERNEL_BASIC)),
		aI({a0,aI}),
		efSPhi(phiS),
		phi(phi),
		nI({n0,nI}),
		j0(j0),
		beta(b),
		n(n)
	{
	}
		
	void FDBVKinetics::init0()
	{
		acl::TypeID type(getElementType(aI[0]->getDContainer()));
	
		acl::VectorOfElements j(acl::generateVEPrivateVariable(1, type));
		auto eta(phi->getEContainer());
		auto phiSL(efSPhi->getEContainer());

		acl::ExpressionContainer kk;
		initMapInfrastructure(kk);

		// valid for beta=0.5 only
		acl::VectorOfElements cProd(aI[0]->getEContainer());
		for(unsigned int i(1); i<aI.size(); ++i)
			copy(cProd*aI[i]->getEContainer(),cProd);
		copy(sqrt(cProd),cProd);
		
		acl::VectorOfElements jCoef(acl::generateVEConstant(0.));
		for(unsigned int i(0); i<aI.size(); ++i)
		{
			acl::VectorOfElements jCoefI(acl::generateVEConstant(1.));
			for(unsigned int j(0); j<aI.size()-1; ++j)
				copy(jCoefI* aI[j+(j>=i?1:0)]->getEContainer() , jCoefI);
			copy(jCoef+ nI[i]*jCoefI,jCoef);
		}
		

		kk << (j = j0*cProd*(exp((1.-beta)*n*eta) - exp(-beta*n*eta))/
		           (1./*-j0/cProd*.5*(exp((1.-beta)*n*eta) - exp(-beta*n*eta))/n*jCoef*/ ));

		//no reaction for 0 anount of reactants
		for(unsigned int i(0); i<aI.size(); ++i)
			if(nI[i]<0)
				kk << (j=min(j,-aI[i]->getEContainer()*(n/nI[i]),type));

		kk << (assignmentSafe(phiSL, j));

		(*kernelJ) << kk;
		(*kernel) << kk;
		for(unsigned int i(0); i<aI.size(); ++i)
			(*kernel) << (assignmentSafe(aI[i]->getEContainer(), 
			                             aI[i]->getEContainer()+nI[i]*j/n));
		
	}

	void FDBVKinetics::addAI(Data ai, double ni)
	{
		aI.push_back(ai);
		nI.push_back(ni);
	}

	void FDBVKinetics::executeJ()
	{
		kernelJ->compute();
	};
		
	SPFDBVKinetics generateFDBVKinetics(SPDataWithGhostNodesACLData a0,
	                                    double n0,
	                                    SPDataWithGhostNodesACLData aI,
	                                    double nI,
	                                    SPDataWithGhostNodesACLData phiS,
	                                    SPAbstractDataWithGhostNodes phi,
	                                    double j0,
	                                    double beta,
	                                    double n)
	{
		auto nm(make_shared<FDBVKinetics> (a0, n0, aI, nI, 
		                                   phiS, phi, 
		                                   acl::generateVEConstant(j0),
		                                   acl::generateVEConstant(beta),
		                                   n));
		return nm;
	}
				
} //asl
