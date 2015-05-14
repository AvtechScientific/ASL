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


#include "aslLBGK.h"
#include "aslGenerators.h"
#include <acl/aclGenerators.h>
#include "acl/acl.h"
#include "acl/DataTypes/aclConstant.h"
#include "math/aslTemplateVE.h"

#include <data/aslDataWithGhostNodes.h>
#include <acl/Kernels/aclKernel.h>
#include <acl/aclElementBase.h>
#include <acl/DataTypes/aclMemBlock.h>
#include <acl/DataTypes/aclIndex.h>
#include <math/aslTemplates.h>
#include "acl/aclMath/aclVectorOfElements.h"
#include "acl/Kernels/aclKernel.h"
#include <math/aslVectorsDynamicLengthOperations.h>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>

namespace asl
{

	acl::VectorOfElements computeRho(acl::VectorOfElements f, const VectorTemplate* vt)
	{
		return f*acl::generateVEConstant(vt->quasiparticlesCoefs);	
	}

	acl::VectorOfElements computeMomentum(acl::VectorOfElements f, const VectorTemplate* vt)
	{
		unsigned int nD(asl::nD(vt->vectors[0]));
		acl::VectorOfElements a(nD);
		unsigned int nv(f.size());
		vector<double> vComponent(nv);
		for(unsigned int i(0); i<nD; ++i){
			for(unsigned int j(0); j<nv; ++j)
				vComponent[j]=double(vt->vectors[j][i])*vt->quasiparticlesCoefs[j];
			a[i]=(f*acl::generateVEConstant(vComponent))[0];
		}
		return a;	
	}

	acl::VectorOfElements a_i_dot_v(acl::VectorOfElements v, const VectorTemplate* vt)
	{
		unsigned int nv(vt->vectors.size());
		acl::VectorOfElements res(nv);
		for(unsigned int j(0); j<nv; ++j)
			res[j]=(v*acl::generateVEConstant(vt->vectors[j]))[0];
		return res;			
	}


	acl::VectorOfElements generateInverceVector(acl::VectorOfElements f, const VectorTemplate* vt)
	{
		unsigned int nv(vt->vectors.size());
		acl::VectorOfElements res(nv);

		for(unsigned int i(0); i<nv; ++i)
			res[i]=f[vt->invertVectors[i]];
		
		return res;
	}

	
	LBGK::LBGK():
		SingleKernelNM(acl::KERNEL_SIMDUA),
		vectorTemplate(NULL),
		flagComputeVelocity(true),
		flagComputeRho(true),
		flagCompressible(true)
	{
	}


	LBGK::LBGK(DataD v, 
	           Param nu, 
	           const VectorTemplate* vT):
		SingleKernelNM(acl::KERNEL_SIMDUA),
		vectorTemplate(vT),
		v(v),
		viscosity(nu),
		flagComputeVelocity(true),
		flagComputeRho(true),
		flagCompressible(true)
	{
		createData(v->getInternalBlock(),
		           v->getEContainer()[0]->getQueue(),
		           v->getEContainer()[0]->getTypeID());
		createCopyKernels();
	}

	LBGK::LBGK(Block b, Param nu, 
		       const VectorTemplate* vT,
	           bool compVel, bool compRho,
		       acl::CommandQueue queue):
		SingleKernelNM(acl::KERNEL_SIMDUA),
		vectorTemplate(vT),
		viscosity(nu),
		flagComputeVelocity(compVel),
		flagComputeRho(compRho),
		flagCompressible(true)
	{
		createData(b,queue,viscosity[0]->getTypeID());
		createCopyKernels();
	}

	void LBGK::createData(Block b, acl::CommandQueue queue, acl::TypeID type)
	{
		unsigned int nv(vectorTemplate->vectors.size());
		fShifts.reset(new AVec<int>(nv));
		fShiftsIncrement.reset(new AVec<int>(nv));
		asl::AVec<int> transformVector(offset(b).c2iTransformVector);
		for(unsigned int i(0); i<nv; ++i){
			(*fShiftsIncrement)[i]=-transformVector*vectorTemplate->vectors[i];
		}

		unsigned int maxIncrement(maxComponent(*fShiftsIncrement));

		unsigned int length(productOfElements(offset(b).getSize()));
		unsigned int poolLength(2*length+4*maxIncrement);

		for(unsigned int i(0); i<nv; ++i){
			(*fShifts)[i]= ((*fShiftsIncrement)[i] >0 ? 0 : poolLength-length);
		}
				
		copy(acl::generateVEData(poolLength, type, nv, queue),fPool);

		acl::VectorOfElements container(generateVESubElements(fPool,
		                                                 length,
		                                                 acl::generateVEVariableSP(fShifts)));
		f=generateDataContainer_SP(b, container, 1u);
		if (!v) 
			v=generateDataContainerACL_SP(b, type, nD(b), 1u, queue);
		if (!rho)
			rho=generateDataContainerACL_SP(b, type, 1u, 1u, queue);
	}

	void LBGK::createCopyKernels()
	{
		unsigned int nv(fPool.size());
		unsigned int length(f->getEContainer()[0]->getSize());
		unsigned int poolLength(fPool[0]->getSize());

		acl::Element ind(new acl::Index(length));
		acl::Element backShift(new acl::Constant<cl_int>(poolLength-length));
		
		copyKernels.resize(nv);
		for(unsigned int i(0); i<nv; ++i){
			copyKernels[i].reset(new acl::Kernel(acl::KERNEL_SIMDUA));
			using namespace acl::elementOperators;
			using acl::elementOperators::operator+;
			if((*fShiftsIncrement)[i]>0)
				copyKernels[i]->addExpression(
				                operatorAssignmentSafe(excerpt(fPool[i], ind),
				                                       f->getEContainer()[i]));
			else
				copyKernels[i]->addExpression(
				                operatorAssignmentSafe(excerpt(fPool[i],backShift+ind),
				                					   f->getEContainer()[i]));

			copyKernels[i]->setup();
		}		
	}

	void LBGK::init0()
	{
		unsigned int nv(vectorTemplate->vectors.size());
		acl::TypeID typeID(viscosity[0]->getTypeID());

		unsigned int nd(v->getEContainer().size());
		
		acl::VectorOfElements fOld(generateVEPrivateVariable(nv,typeID));
		acl::VectorOfElements fEq(generateVEPrivateVariable(nv,typeID));
		acl::VectorOfElements pRho(generateVEPrivateVariable(1u, typeID));

		acl::VectorOfElements momentum(generateVEPrivateVariable(nd,typeID));
		acl::VectorOfElements aDotMomentum(a_i_dot_v(momentum,vectorTemplate));
		
		acl::VectorOfElements w(generateVEPrivateVariable(1u, typeID));
		(*kernel)<<(fOld=f->getEContainer());
		(*kernel)<<(w= 1./(3. * viscosity + .5));
		(*kernel)<<(momentum = computeMomentum(fOld,vectorTemplate));
		(*kernel)<<(pRho = computeRho(fOld,vectorTemplate));
		// valid for both compressible and incompresible
		if(omega.size()>0)
			(*kernel)<< (momentum += (2./w)* crossProduct(omega, momentum));

		if(flagComputeRho)
			(*kernel)<<	acl::assignmentSafe(rho->getEContainer(), pRho);
		if(flagComputeVelocity)
		{
			if (flagCompressible)
				(*kernel)<<	acl::assignmentSafe(v->getEContainer(), momentum/pRho);
			else
				(*kernel)<<	acl::assignmentSafe(v->getEContainer(), momentum);
		}
		
		if (flagCompressible)
			(*kernel)<<(fEq=pRho * acl::generateVEConstantN(nv, 1.) + 
			                3. * aDotMomentum +
				            4.5 * acl::productOfElements(aDotMomentum,aDotMomentum)/pRho -
		    		    	1.5 * momentum*momentum/pRho * acl::generateVEConstantN(nv, 1.)
		        		);
		else
			(*kernel)<<(fEq=pRho * acl::generateVEConstantN(nv, 1.) + 
			                3. * aDotMomentum +
				            4.5 * acl::productOfElements(aDotMomentum,aDotMomentum) -
		    		    	1.5 * momentum*momentum * acl::generateVEConstantN(nv, 1.)
		        		);
			
		auto fNew(fOld*(1.-w)+ w*fEq);
		(*kernel)<<(assignmentSafe(f->getEContainer(),fNew));		
	}

	void LBGK:: preProcessing()
	{
		unsigned int nv(fShifts->getSize());
		unsigned int length(f->getEContainer()[0]->getSize());
		unsigned int poolLength(fPool[0]->getSize());
		for(unsigned int i(0); i<nv; ++i)
			if((*fShifts)[i]+length+(*fShiftsIncrement)[i]>poolLength || 
			   (*fShifts)[i]+(*fShiftsIncrement)[i]<0)
			{
				copyKernels[i]->compute();
				(*fShifts)[i]=(*fShiftsIncrement)[i] >0 ? 0 : poolLength-length;
			}
		(*fShifts)+=(*fShiftsIncrement);
	}

	void LBGK::setOmega(LBGK::Param w)
	{
		copy(w,omega);
	}

	LBGKUtilities::LBGKUtilities(SPLBGK lbgk):
		num(lbgk)
	{
	}

	void LBGKUtilities::initF(Param rho, Param vel)
	{
		unsigned int nv(num->getVectorTemplate()->vectors.size());
		acl::VectorOfElements aDotVel(a_i_dot_v(vel,num->getVectorTemplate()));

		initData(num->getF()->getEContainer(),
		         rho* (acl::generateVEConstantN(nv, 1.) + 
		                3. * aDotVel +
		                4.5 * acl::productOfElements(aDotVel,aDotVel) -
		            	1.5 * vel * vel * acl::generateVEConstantN(nv, 1.)),
		         acl::KERNEL_SIMDUA);
		initData(num->getRho()->getEContainer(),rho);
		initData(num->getVelocity()->getEContainer(),vel);
	}

	void LBGKUtilities::initF(Param vel)
	{
		initF(acl::generateVEConstant(1),vel);
	}	

	LBGKTurbulence::LBGKTurbulence(Block b, Param nu, 
	                               const VectorTemplate* vT,
	                               bool compVel, bool compRho,
	                               acl::CommandQueue queue):
		LBGK(b,nu,vT,compVel, compRho, queue)
	{
	}

	LBGKTurbulence::LBGKTurbulence(DataD v, 
	                               Param nu, 
	                               const VectorTemplate* vT):
		LBGK(v,nu,vT)
	{
	}


	void LBGKTurbulence::init0()
	{
		unsigned int nv(vectorTemplate->vectors.size());
		acl::TypeID typeID(viscosity[0]->getTypeID());

		unsigned int nd(v->getEContainer().size());
		
		acl::VectorOfElements fOld(generateVEPrivateVariable(nv,typeID));
		acl::VectorOfElements fEq(generateVEPrivateVariable(nv,typeID));
		acl::VectorOfElements pRho(generateVEPrivateVariable(1u, typeID));

		acl::VectorOfElements momentum(generateVEPrivateVariable(nd,typeID));
		acl::VectorOfElements aDotMomentum(a_i_dot_v(momentum,vectorTemplate));
		
		acl::VectorOfElements w(generateVEPrivateVariable(1u, typeID));
		(*kernel)<<(fOld=f->getEContainer());
// Turbulence
		auto fL2(fOld*fOld);
		auto fmfEqL2((fOld-fEq)*(fOld-fEq));
		(*kernel)<<(w=4.*(fmfEqL2/fL2));
		(*kernel)<<(w = select((1./(3. * viscosity + .5)),
		                       (1./(3. * min(viscosity*w,
		                                     acl::generateVEConstant(.5)) + .5)), 
		                       w>1,
		                       typeID));
//--------------------------
//		(*kernel)<<(w= 1./(3. * viscosity + .5));
		(*kernel)<<(momentum = computeMomentum(fOld,vectorTemplate));
		(*kernel)<<(pRho = computeRho(fOld,vectorTemplate));
		// valid for both compressible and incompresible
		if(omega.size()>0)
			(*kernel)<< (momentum += (2./w)* crossProduct(omega, momentum));

		if(flagComputeRho)
			(*kernel)<<	acl::assignmentSafe(rho->getEContainer(), pRho);
		if(flagComputeVelocity)
		{
			if (flagCompressible)
				(*kernel)<<	acl::assignmentSafe(v->getEContainer(), momentum/pRho);
			else
				(*kernel)<<	acl::assignmentSafe(v->getEContainer(), momentum);
		}
		
		if (flagCompressible)
			(*kernel)<<(fEq=pRho * acl::generateVEConstantN(nv, 1.) + 
			                3. * aDotMomentum +
				            4.5 * acl::productOfElements(aDotMomentum,aDotMomentum)/pRho -
		    		    	1.5 * momentum*momentum/pRho * acl::generateVEConstantN(nv, 1.)
		        		);
		else
			(*kernel)<<(fEq=pRho * acl::generateVEConstantN(nv, 1.) + 
			                3. * aDotMomentum +
				            4.5 * acl::productOfElements(aDotMomentum,aDotMomentum) -
		    		    	1.5 * momentum*momentum * acl::generateVEConstantN(nv, 1.)
		        		);
		auto fNew(fOld*(1.-w)+ w*fEq);
		(*kernel)<<(assignmentSafe(f->getEContainer(),fNew));		
	}
		
} //asl

