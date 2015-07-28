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
	\example testVectorOfElements.cc
 */

#include "acl/acl.h"
#include "utilities/aslUValue.h"
#include "acl/aclGenerators.h"
#include "acl/Kernels/aclKernel.h"
#include <acl/aclMath/aclVectorOfElements.h>
#include "aslUtilities.h"
#include <math.h>

using namespace asl;
using namespace acl;
using namespace std;

bool testSimpleKernel()
{
	cout << "Test of \"Simple kernel\" function..." << flush;

	VectorOfElements vec0(3);
	VectorOfElements vec1(1);
	copy(generateVEData<float>(10u,3u),vec0);
	copy(generateVEData<float>(10u,1u),vec1);

	Kernel k;
	{ 
		using namespace elementOperators;
		k << (vec0=generateVEConstant(0.1f,1.f,2.f));
		k << (vec1=generateVEConstant(0.f,1.f,2.f)*vec0+generateVEConstant(1.f));
	}
	k.setup();
	k.compute();

	vector<cl_float> output0(10), output1(10), output2(10);
	vector<cl_float> output3(10);
	copy(vec0[0], output0);	
	copy(vec0[1], output1);
	copy(vec0[2], output2);
	copy(vec1[0], output3);
	
	bool status(output0[9]<0.101 && output1[2] ==1 && output2[5] ==2. && output3[1] ==6.);
	errorMessage(status);

	return status;		
}


bool testAdvancedOperations()
{
	cout << "Test of advanced operations..." << flush;
	VectorOfElements vec0(2);
	VectorOfElements vec1(2);
	VectorOfElements res(2);
	VectorOfElements ind(1);
	VectorOfElements c(1);
	
	asl::UValue<cl_int> v0(7);
	asl::UValue<cl_int> v1(8);
	copy(generateVEConstant(3, 4), vec0);
	copy(generateVEVariableSP(v0.p, v1.p), vec1);
	copy(generateVEData<int>(11, 2), res);
	copy(generateVEIndex(), ind);
	copy(generateVEConstant(5), c);

	vector<cl_int> output0(11, 0);
	vector<cl_int> output1(11, 3);
	vector<cl_int> expected0({10, 10, 10, 10, 10, 10, -4, -4, -4, -4, -4});
	vector<cl_int> expected1({12, 12, 12, 12, 12, 12, -4, -4, -4, -4, -4});


	Kernel k;
	{ 
		k << (res = select(vec0 + vec1, vec0 - vec1, ind > c));
	}
	k.setup();
	k.compute();

	copy(res[0], output0);
	copy(res[1], output1);
	
	bool status(output0 == expected0 && output1 == expected1);
	errorMessage(status);

	return status;
}


int main()
{
	bool allTestsPassed(true);

	allTestsPassed &= testSimpleKernel();
	allTestsPassed &= testAdvancedOperations();

	return allTestsPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}