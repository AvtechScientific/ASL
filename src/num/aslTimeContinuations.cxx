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


#include "aslTimeContinuations.h"
#include <aslGenerators.h>
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include "aslUtilities.h"
#include <acl/aclMath/aclVectorOfElements.h>

#include <data/aslDataWithGhostNodes.h>
#include <acl/Kernels/aclKernel.h>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>

namespace asl
{

	TimeContinuations::TimeContinuations(Data inD, double f):
		inData(inD->getDContainer()),
		factor(f),
		nStorages(0)		
	{
	}

	TimeContinuations::TimeContinuations(acl::VectorOfElementsData & inD, double f):
		inData(inD),
		factor(f),
		nStorages(0)		
	{
	}

	void TimeContinuations::addData(acl::VectorOfElementsData & inD)
	{
		copy(cat(inData,inD),inData);
	}	

	void TimeContinuations::addData(Data inD)
	{
		copy(cat(inData,inD->getDContainer()),inData);
	}	

	void TimeContinuations::reset()
	{
		nStorages=0;
	}			
		
	double polynomRootsLag(double x, int pOrder, int iCoef)
	{
		double res(1);
		for(int i(0); i<pOrder; ++i)
		{
			int a= (i>=iCoef) ? 1 : 0;
			res*=x+(pOrder-i-a);
		}					
		return res; 
	}
	
	double polynomCoefsLag(double f, int pOrder, int iCoef)
	{
		return polynomRootsLag(f,pOrder,iCoef)/polynomRootsLag(iCoef-pOrder,pOrder,iCoef);		
	}
	

	TimeContinPLagrange::TimeContinPLagrange(Data inD, double f, unsigned int o):
		TimeContinuations(inD,f),
		kernels(o+1),
		order(o),
		coefs(o+1)
	{
	}

	TimeContinPLagrange::TimeContinPLagrange(acl::VectorOfElementsData & inD, 
	                                         double f, unsigned int o):
		TimeContinuations(inD,f),
		kernels(o+1),
		order(o),
		coefs(o+1)
	{
	}
		
	void TimeContinPLagrange::init()
	{
		unsigned int nC(inData.size());
		copy(clone(inData),storedData);
		acl::initData(storedData, acl::generateVEConstantN(nC,0.));

		acl::VectorOfElements inDataV(inData);

		for(unsigned int i(0); i<=order; ++i)
			coefs[i]=polynomCoefsLag(factor, order, i);

//		cout<<factor<<endl;
//		cout<<coefs<<endl;

		for(unsigned int i(0); i<=order; ++i)
			kernels[i].reset(new acl::Kernel(acl::KERNEL_SIMD));

		(*kernels[0]) << (assignmentSafe(storedData,coefs[0]*inDataV));
		for(unsigned int i(1); i<order; ++i)
			(*kernels[i]) << (assignmentSafe(storedData,storedData+coefs[i]*inDataV));
		
		(*kernels[order])<< (assignmentSafe (inDataV, storedData+coefs[order]*inDataV));

		setupAll(kernels);
	}
			
	void TimeContinPLagrange::execute()
	{
		kernels[nStorages]->compute();
		nStorages++;
		nStorages=nStorages % (order+1);
	}	

	double polynomRootsLagFr(double x, double offset, int pOrder, int iCoef)
	{
		double res(1);
		for(int i(0); i<pOrder; ++i)
		{
			int a= (i>=iCoef) ? 1 : 0;
			res*=1./(offset+x)-1./(offset-pOrder+i+a);
		}					
		return res; 
	}
	
	double polynomCoefsLagFr(double f,double offset, int pOrder, int iCoef)
	{
		return polynomRootsLagFr(f,offset,pOrder,iCoef)/
			   polynomRootsLagFr(iCoef-pOrder,offset,pOrder,iCoef);		
	}
	

	TimeContinPLagrangeFraction::TimeContinPLagrangeFraction(Data inD, double f, unsigned int o):
		TimeContinuations(inD,f),
		kernels(o+1),
		order(o),
//		offset(.5*factor+o*2),
		offset(1.+o*2),
		coefs(o+1)
	{
	}

	TimeContinPLagrangeFraction::TimeContinPLagrangeFraction(acl::VectorOfElementsData & inD, 
	                                                         double f, unsigned int o):
		TimeContinuations(inD,f),
		kernels(o+1),
		order(o),
//		offset(.5*factor+o*2),
		offset(1.+o*2),
		coefs(o+1)
	{
	}
		
	void TimeContinPLagrangeFraction::init()
	{
		unsigned int nC(inData.size());
		copy(clone(inData),storedData);
		acl::initData(storedData, acl::generateVEConstantN(nC,0.));
		acl::VectorOfElements inDataV(inData);

		for(unsigned int i(0); i<=order; ++i)
			coefs[i]=polynomCoefsLagFr(factor, offset, order, i);

//		cout<<factor<<endl;
//		cout<<coefs<<endl;

		for(unsigned int i(0); i<=order; ++i)
			kernels[i].reset(new acl::Kernel(acl::KERNEL_SIMD));

		(*kernels[0]) << (assignmentSafe(storedData,coefs[0]*inDataV));
		for(unsigned int i(1); i<order; ++i)
			(*kernels[i]) << (assignmentSafe(storedData,storedData+coefs[i]*inDataV));
		
		(*kernels[order])<< (assignmentSafe (inDataV, storedData+coefs[order]*inDataV));

		setupAll(kernels);
	}
			
	void TimeContinPLagrangeFraction::execute()
	{
		kernels[nStorages]->compute();
		nStorages++;
		nStorages=nStorages % (order+1);
	}	
		
} //asl
