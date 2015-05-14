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


#include "aslFDPoroElasticity.h"
#include <aslGenerators.h>
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <acl/Kernels/aclKernel.h>
#include <math/aslTemplateVE.h>
#include <math/aslTemplates.h>
#include <data/aslDataWithGhostNodes.h>
#include <aslDataInc.h>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>

namespace asl
{
	FDPoroElasticity::FDPoroElasticity():
		ElasticityCommonA(),
		pressureData(),
		pressureInternalData(),
		pressureLiquidData(),
		pressureLiquidInternalData(),		
		nSubsteps(acl::generateVEConstant(10.))
	{
	}

 	FDPoroElasticity::~FDPoroElasticity()
	{
	}
		
	FDPoroElasticity::FDPoroElasticity(Data d, Data pl, Param bM, 
	                                   Param sM, Param k,
	                                   const VectorTemplate* vT):
		ElasticityCommonA(d,bM,sM,vT),		
		pressureData(),
		pressureInternalData(),
		pressureLiquidData(pl),
		pressureLiquidInternalData(),
		hydraulicCondactivity(k),
		nSubsteps(acl::generateVEConstant(10.))
	{		
	}

	void FDPoroElasticity::setNSubsteps(unsigned int n)
	{
		nSubsteps=acl::generateVEConstant(n);
	}
	
	void FDPoroElasticity::init()
	{

		unsigned int nD(displacementData->getDContainer().size());
		displacementInternalData=clone(displacementData);
		pressureData = clone(pressureLiquidData);
		pressureInternalData = clone(pressureLiquidData);
		pressureLiquidInternalData = clone(pressureLiquidData);

		initData(displacementInternalData->getEContainer(), 
		         displacementData->getEContainer());
		asl::initData(pressureData, 0.);
		asl::initData(pressureInternalData, 0.);
		initData(pressureLiquidInternalData->getEContainer(), 
		         pressureLiquidData->getEContainer());

		acl::TypeID type(getElementType(displacementData->getDContainer()));
		         
		acl::VectorOfElements dnew(acl::generateVEPrivateVariable(3,type)); 
		acl::VectorOfElements pnew(acl::generateVEPrivateVariable(1,type));
		acl::VectorOfElements plnew(acl::generateVEPrivateVariable(1,type));
		acl::VectorOfElements coef1(acl::generateVEPrivateVariable(1,type));
		acl::VectorOfElements coef2(acl::generateVEPrivateVariable(1,type));

		double sFactor(10.);         //< artificial compresibility
		acl::VectorOfElements coefPlP(acl::generateVEConstant(sFactor));

		static const double w(.5); 
		
		acl::VectorOfElements wC(acl::generateVEConstant(w));
		acl::VectorOfElements wm1C(acl::generateVEConstant(1.-w));
		
		(*kernel)<<(coef1 = shearModulus);
		(*kernel)<<(coef2 = (bulkModulus+shearModulus/3.));

		TemplateVE du;
		du.init(*pressureData,*vectorTemplate);
		(*kernel)<< du.initValues << (pnew=du.getValue(0)*wm1C);
		(*kernel)<<(dnew=-gradient(du));
		(*kernel)<<(plnew=-coefPlP*du.getValue(0));

		du.init(*pressureLiquidData,*vectorTemplate);
		(*kernel)<< du.initValues; 
		(*kernel)<< (plnew+=du.getValue(0) + hydraulicCondactivity/nSubsteps*laplas(du));
		(*kernel)<<(dnew-=gradient(du));
				
		for (unsigned int i(0); i<nD; ++i)
		{	
			du.init(*displacementData,*vectorTemplate,i);
			(*kernel)<< du.initValues;
			(*kernel)<< (pnew -= coef2*subVE(gradient(du),i)*wC);
			(*kernel)<< (subVE(dnew,i) += du.getValue(0) + coef1 * laplas(du));
		}
		(*kernel)<<(plnew+=coefPlP*pnew);

		if(force.size()>0)
			(*kernel)<<(dnew+= force);

		(*kernel)<<(assignmentSafe (displacementInternalData->getSubContainer(),dnew));
		(*kernel)<<(assignmentSafe (pressureInternalData->getSubContainer(),pnew));
		(*kernel)<<(assignmentSafe (pressureLiquidInternalData->getSubContainer(),plnew));

		kernel->setup();

	}


	void FDPoroElasticity::execute()
	{
		kernel->compute();
		swapBuffers(displacementData->getDContainer(),
		            displacementInternalData->getDContainer());
		swapBuffers(pressureData->getDContainer(), 
		            pressureInternalData->getDContainer());
		swapBuffers(pressureLiquidData->getDContainer(), 
		            pressureLiquidInternalData->getDContainer());
	}

		SPFDPoroElasticity generateFDPoroElasticity(SPDataWithGhostNodesACLData d,
		                                            SPDataWithGhostNodesACLData pl,
		                                            double bM, 
		                                            double sM, 
		                                            double k,
		                                            const VectorTemplate* vT)
	{
		auto nm(make_shared<FDPoroElasticity> (d, 
		                                       pl,
		                                       acl::generateVEConstant(bM),
		                                       acl::generateVEConstant(sM),
		                                       acl::generateVEConstant(k),		                                       
		                                       vT));
		return nm;				
	}
		
	template <typename T> 
		SPFDPoroElasticity generateFDPoroElasticity(SPDataWithGhostNodesACLData d,
		                                            SPDataWithGhostNodesACLData pl,
		                                            UValue<T> bM, 
		                                            UValue<T> sM,
		                                            UValue<T> k,
		                                            const VectorTemplate* vT)
	{
		auto nm(make_shared<FDPoroElasticity> (d, 
		                                       pl,
		                                       acl::generateVEVariableSP(bM.p),
		                                       acl::generateVEVariableSP(sM.p),
		                                       acl::generateVEVariableSP(k.p),		                                       
		                                       vT));
		return nm;						
	}

	template  
		SPFDPoroElasticity generateFDPoroElasticity(SPDataWithGhostNodesACLData d,
		                                            SPDataWithGhostNodesACLData pl,
		                                            UValue<double> bM, 
		                                            UValue<double> sM,
		                                            UValue<double> k,
		                                            const VectorTemplate* vT);
	template  
		SPFDPoroElasticity generateFDPoroElasticity(SPDataWithGhostNodesACLData d,
		                                            SPDataWithGhostNodesACLData pl,
		                                            UValue<float> bM, 
		                                            UValue<float> sM,
		                                            UValue<float> k,
		                                            const VectorTemplate* vT);			
		
} //asl
