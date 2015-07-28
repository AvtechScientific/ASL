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


/**
 	\example testPrivateVar.cc
*/

#include "utilities/aslUValue.h"
#include "acl/Kernels/aclKernel.h"
#include "acl/Kernels/aclKernelConfigurationTemplates.h"
#include <acl/acl.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <acl/aclGenerators.h>
#include <utilities/aslTimer.h>

const unsigned int nLength(1000000);
const unsigned int nOperations(10);
const unsigned int nCycles(10000);
const acl::KernelConfiguration & kConf(acl::KERNEL_BASIC); 
//const acl::KernelConfiguration & kConf(acl::KERNEL_SIMDUA); 

bool testKernelUnoptimized()
{
	cout << "Test of \"Simple kernel\" function..." << flush;

	auto vec1(acl::generateVEData<float>(nLength,1u));
	acl::Kernel k(kConf);

	auto res(acl::generateVEPrivateVariable(1u, acl::TYPE_FLOAT));
	k << (res = acl::generateVEConstant(0.));
	for(unsigned int i(0); i<nOperations; ++i)
	{
		auto tempRes(acl::generateVEPrivateVariable(1u, acl::TYPE_FLOAT));
		k << (tempRes=acl::generateVEIndex() * i+5.*i*i+7./(i+1.));
		k << (res+=tempRes*tempRes+rsqrt(.2*tempRes)+ 1./tempRes);
	}
	k << (vec1=res);
	k.setup();

	asl::Timer timer;
	timer.start();
	for(unsigned int i(0); i<nCycles; ++i)
		k.compute();
	timer.stop();
	std::cout<<"Unoptimized: "<<timer.getTime()<<endl;

	return true;		
}


bool testKernelUnoptimizedPlus()
{
	cout << "Test of \"Simple kernel\" function..." << flush;

	auto vec1(acl::generateVEData<float>(nLength,1u));
	acl::Kernel k(kConf);

	auto res(acl::generateVEPrivateVariable(1u, acl::TYPE_FLOAT));
	k << (res = acl::generateVEConstant(0.));
	vector<acl::VectorOfElements> tempRes(nOperations);
	for(unsigned int i(0); i<nOperations; ++i)
	{
		copy(acl::generateVEPrivateVariable(1u, acl::TYPE_FLOAT), tempRes[i]);
		k << (tempRes[i]=acl::generateVEIndex() * i+5.*i*i+7./(i+1.));
	}
	for(unsigned int i(0); i<nOperations; ++i)
		k << (res+=tempRes[i]*tempRes[i]+rsqrt(.2*tempRes[i])+ 1./tempRes[i]);
	k << (vec1=res);
	k.setup();

	asl::Timer timer;
	timer.start();
	for(unsigned int i(0); i<nCycles; ++i)
		k.compute();
	timer.stop();
	std::cout<<"UnoptimizedPlus: "<<timer.getTime()<<endl;

	return true;		
}


bool testKernelOptimized()
{
	cout << "Test of \"Simple kernel\" function..." << flush;

	auto vec1(acl::generateVEData<float>(nLength,1u));
	acl::Kernel k(kConf);

	auto res(acl::generateVEPrivateVariable(1u, acl::TYPE_FLOAT));
	auto tempRes(acl::generateVEPrivateVariable(1u, acl::TYPE_FLOAT));
	k << (res = acl::generateVEConstant(0.));
	for(unsigned int i(0); i<nOperations; ++i)
	{
		k << (tempRes=acl::generateVEIndex() * i+5.*i*i+7./(i+1.));
		k << (res+=tempRes*tempRes+rsqrt(.2*tempRes)+ 1./tempRes);
	}
	k << (vec1=res);
	k.setup();
//	cout<<k.getKernelSource()<<endl;

	asl::Timer timer;
	timer.start();
	for(unsigned int i(0); i<nCycles; ++i)
		k.compute();
	timer.stop();
	std::cout<<"Optimized: "<<timer.getTime()<<endl;

	return true;		
}



int main()
{
	testKernelUnoptimized();
	testKernelUnoptimizedPlus();
	testKernelOptimized();

	return 0;
}
