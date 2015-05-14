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


#include "aslFDElasticity.h"
#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <acl/Kernels/aclKernel.h>
#include <math/aslTemplateVE.h>
#include <math/aslTemplates.h>
#include <data/aslDataWithGhostNodes.h>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>

namespace asl
{

	ElasticityCommonA::ElasticityCommonA():
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		displacementData(),
		displacementInternalData(),
		vectorTemplate(NULL)
	{
	}

 	ElasticityCommonA::~ElasticityCommonA()
	{
	}
		
	ElasticityCommonA::ElasticityCommonA(Data d, 
	                                     Param bM, 
	                                     Param sM, 
	                                     const VectorTemplate* vT):
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		displacementData(d),
		displacementInternalData(),
		bulkModulus(bM),
		shearModulus(sM),
		vectorTemplate(vT)
	{		
	}

	void ElasticityCommonA::setForce(FDElasticityIncompressibleStatic::Param f)
	{
		acl::copy(f,force);
	}

	unsigned int nD(const ElasticityCommonA & e)
	{
		return nD(e.getDisplacementData()->getBlock());
	}
		
	FDElasticityIncompressibleStatic::FDElasticityIncompressibleStatic():
		ElasticityCommonA(),
		pressure(),
		pressureInternalData()		
	{
	}

 	FDElasticityIncompressibleStatic::~FDElasticityIncompressibleStatic()
	{
	}
		
	FDElasticityIncompressibleStatic::FDElasticityIncompressibleStatic(Data d, 
	                                                       Param bM, 
	                                                       Param sM, 
	                                                       const VectorTemplate* vT):
		ElasticityCommonA(d,bM,sM,vT),
		pressure(),
		pressureInternalData()
	{		
	}

	
	void FDElasticityIncompressibleStatic::init(){
		unsigned int nD(displacementData->getDContainer().size());
		displacementInternalData = clone(displacementData);
		pressure = clone(displacementData, 1u);
		pressureInternalData = clone(pressure);
		acl::initData(displacementInternalData->getDContainer(),
		              displacementData->getEContainer());
		asl::initData(pressure, 0.);
		asl::initData(pressureInternalData, 0.);

		acl::TypeID type(getElementType(displacementData->getDContainer()));
		         
		acl::VectorOfElements dnew(acl::generateVEPrivateVariable(3,type)); 
		acl::VectorOfElements pnew(acl::generateVEPrivateVariable(1,type)); 
		acl::VectorOfElements coef1(acl::generateVEPrivateVariable(1,type));
		acl::VectorOfElements coef2(acl::generateVEPrivateVariable(1,type));

		static const double w(.5); 
		acl::VectorOfElements wC(acl::generateVEConstant(w));
		acl::VectorOfElements wm1C(acl::generateVEConstant(1.-w));
		(*kernel)<<(coef1 = shearModulus);
		(*kernel)<<(coef2 = (bulkModulus+shearModulus/3.));

		TemplateVE du;
		du.init(*pressure,*vectorTemplate);
		(*kernel)<< du.initValues << (pnew=du.getValue(0)*wm1C);
		(*kernel)<<(dnew=displacementData->getSubContainer()-gradient(du));
		
		auto lap(acl::generateVEPrivateVariable(nD,type));
		
		for (unsigned int i(0); i<nD; ++i)
		{	
			du.init(*displacementData,*vectorTemplate,i);
			(*kernel)<< du.initValues << (subVE(lap,i) = laplas(du));
			(*kernel)<< (pnew -= coef2*subVE(gradient(du),i)*wC);
			(*kernel)<< (subVE(dnew,i) += coef1 * subVE(lap,i));
//			for (unsigned int k(0); k<nD; ++k)
//				(*kernel)<<(subVE(dnew,k) += coef1 * dIdJ(k,i,du)); //this commment improoves efficiency
//			kernel->addExpression(acl::elementOperators::barrier());
		}

		if(force.size()>0)
			(*kernel)<<(dnew+= force);

		(*kernel)<<(assignmentSafe (displacementInternalData->getSubContainer(),dnew));
		(*kernel)<<(assignmentSafe (pressureInternalData->getSubContainer(),pnew));

		kernel->setup();
	}


	void FDElasticityIncompressibleStatic::execute()
	{
		kernel->compute();
		swapBuffers(displacementData->getDContainer(),
		            displacementInternalData->getDContainer());
		swapBuffers(pressure->getDContainer(),
		            pressureInternalData->getDContainer());
	}

	FDElasticityRelaxation::FDElasticityRelaxation():
		ElasticityCommonA(),
		pressure(),
		pressureInternalData()		
	{
	}

 	FDElasticityRelaxation::~FDElasticityRelaxation()
	{
	}
		
	FDElasticityRelaxation::FDElasticityRelaxation(Data d, 
	                                               Param bM, 
	                                               Param sM,
	                                               Param dt,
	                                               const VectorTemplate* vT):
		ElasticityCommonA(d,bM,sM,vT),
		pressure(),
		pressureInternalData(),
		deltat(dt)
	{		
	}

	
	void FDElasticityRelaxation::init(){
		unsigned int nD(displacementData->getDContainer().size());
		displacementInternalData = clone(displacementData);
		pressure = clone(displacementData, 1u);
		pressureInternalData = clone(pressure);
		acl::initData(displacementInternalData->getDContainer(),
		              displacementData->getEContainer());
		asl::initData(pressure, 0.);
		asl::initData(pressureInternalData, 0.);

		acl::TypeID type(getElementType(displacementData->getDContainer()));
		         
		acl::VectorOfElements dnew(acl::generateVEPrivateVariable(3,type)); 
		acl::VectorOfElements pnew(acl::generateVEPrivateVariable(1,type)); 
		acl::VectorOfElements coef1(acl::generateVEPrivateVariable(1,type));
		acl::VectorOfElements coef2(acl::generateVEPrivateVariable(1,type));

		static const double w(.5); 
		acl::VectorOfElements wC(acl::generateVEConstant(w));
		acl::VectorOfElements wm1C(acl::generateVEConstant(1.-w));
		(*kernel)<<(coef1=deltat*deltat*shearModulus);
		(*kernel)<<(coef2=(deltat*deltat)*(bulkModulus+shearModulus/3.));

		TemplateVE du;
		du.init(*pressure,*vectorTemplate);
		(*kernel)<< du.initValues << (pnew=du.getValue(0)*wm1C);
		if(dumpingFactor.size()==0)
			(*kernel)<<(dnew=2. * displacementData->getSubContainer()
			                 - displacementInternalData->getSubContainer());
		else
			(*kernel)<<(dnew=(1.+dumpingFactor) * displacementData->getSubContainer()
			                 -dumpingFactor*displacementInternalData->getSubContainer());
//		(*kernel)<<(dnew+=.2* gradient(du)); //!!! should be -
		
		auto lap(acl::generateVEPrivateVariable(nD,type));
		
		for (unsigned int i(0); i<nD; ++i)
		{	
			du.init(*displacementData,*vectorTemplate,i);
			(*kernel)<< du.initValues << (subVE(lap,i) = laplas(du));
			(*kernel)<< (pnew -= coef2*subVE(gradient(du),i)*wC);
			(*kernel)<< (subVE(dnew,i) += coef1 * subVE(lap,i));
			for (unsigned int k(0); k<nD; ++k)
				(*kernel)<<(subVE(dnew,k) += coef2 * dIdJ(k,i,du)); //!!!should be commentent
			kernel->addExpression(acl::elementOperators::barrier());
		}

		if(force.size()>0)
			(*kernel)<<(dnew+= force);

		(*kernel)<<(assignmentSafe (displacementInternalData->getSubContainer(),dnew));
		(*kernel)<<(assignmentSafe (pressureInternalData->getSubContainer(),pnew));

		kernel->setup();
	}


	void FDElasticityRelaxation::execute()
	{
		kernel->compute();
		swapBuffers(displacementData->getDContainer(),
		            displacementInternalData->getDContainer());
		swapBuffers(pressure->getDContainer(),
		            pressureInternalData->getDContainer());
	}

	void FDElasticityRelaxation::setDumpingFactor(FDElasticity2::Param dumpF)
	{
		acl::copy(dumpF,dumpingFactor);
	}

		
	FDElasticity2::FDElasticity2():
		ElasticityCommonA()
//		velocityData(),
	{
	}

	FDElasticity2::FDElasticity2(Data d, Param bM, 
	                             Param sM, Param dt, 
	                             const VectorTemplate* vT):
		ElasticityCommonA(d,bM,sM,vT),
		deltat(dt)
	{		
	}

 	FDElasticity2::~FDElasticity2()
	{
	}
		
	void FDElasticity2::init()
	{
		unsigned int nD(displacementData->getDContainer().size());
		displacementInternalData=clone(displacementData);
		acl::initData(displacementInternalData->getDContainer(),
		              displacementData->getEContainer());

		acl::TypeID type(getElementType(displacementData->getDContainer()));
		         
		acl::VectorOfElements dnew(acl::generateVEPrivateVariable(3,type)); 
		acl::VectorOfElements coef1(acl::generateVEPrivateVariable(1,type));
		acl::VectorOfElements coef2(acl::generateVEPrivateVariable(1,type));
		(*kernel)<<(coef1=deltat*deltat*shearModulus);
		(*kernel)<<(coef2=(deltat*deltat)*(bulkModulus+shearModulus/3.));
		if(dumpingFactor.size()==0)
			(*kernel)<<(dnew=2. * displacementData->getSubContainer()
			                 - displacementInternalData->getSubContainer());
		else
			(*kernel)<<(dnew=(1.+dumpingFactor) * displacementData->getSubContainer()
			                 -dumpingFactor*displacementInternalData->getSubContainer()
		               );

		TemplateVE du;
		auto lap(acl::generateVEPrivateVariable(nD,type));

		for (unsigned int i(0); i<nD; ++i)
		{	
			du.init(*displacementData,*vectorTemplate,i);
			(*kernel)<< du.initValues << (subVE(lap,i) = laplas(du));
			(*kernel)<< (subVE(dnew,i) += coef1 * subVE(lap,i));
			for (unsigned int k(0); k<nD; ++k)
				(*kernel)<<(subVE(dnew,k) += coef2 * dIdJ(k,i,du));
//			kernel->addExpression(acl::elementOperators::barrier());
		}

		if(force.size()>0)
			(*kernel)<<(dnew+= deltat * deltat * force);

		(*kernel)<<(assignmentSafe (displacementInternalData->getSubContainer(),dnew));

		kernel->setup();
	}


	void FDElasticity2::execute()
	{
		kernel->compute();
		swapBuffers(displacementData->getDContainer(),displacementInternalData->getDContainer());
	}
	
	void FDElasticity2::setDumpingFactor(FDElasticity2::Param dumpF)
	{
		acl::copy(dumpF,dumpingFactor);
	}

	SPFDElasticity2 generateFDElasticity(SPDataWithGhostNodesACLData d, 
	                                     double bM, 
	                                     double sM, 
	                                     double dt, 
	                                     const VectorTemplate* vT)
	{
		auto nm(make_shared<FDElasticity2> (d, 
		                                    acl::generateVEConstant(bM),
		                                    acl::generateVEConstant(sM),
		                                    acl::generateVEConstant(dt),
		                                    vT));
		return nm;		
	}

	SPFDElasticityIncompressibleStatic generateFDElasticityStatic(SPDataWithGhostNodesACLData d, 
	                                                              double bM, 
	                                                              double sM, 
	                                                              const VectorTemplate* vT)
	{
		auto nm(make_shared<FDElasticityIncompressibleStatic> (d, 
		                                                       acl::generateVEConstant(bM),
		                                                       acl::generateVEConstant(sM),
		                                                       vT));
		return nm;				
	}

	SPFDElasticityRelaxation generateFDElasticityRelax(SPDataWithGhostNodesACLData d, 
	                                                   double bM, 
	                                                   double sM, 
	                                                   double dt, 
	                                                   const VectorTemplate* vT)
	{
		auto nm(make_shared<FDElasticityRelaxation> (d, 
		                                             acl::generateVEConstant(bM),
		                                             acl::generateVEConstant(sM),
		                                             acl::generateVEConstant(dt),
		                                             vT));
		return nm;		
	}
		
} //asl
