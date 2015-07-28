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
	\example testKernel.cc
 */


#include "acl/acl.h"
#include "acl/DataTypes/aclIndex.h"
#include "acl/DataTypes/aclGroupID.h"
#include "acl/DataTypes/aclConstant.h"
#include "acl/DataTypes/aclVariable.h"
#include "acl/DataTypes/aclVariableReference.h"
#include "acl/DataTypes/aclPrivateVariable.h"
#include "acl/DataTypes/aclPrivateArray.h"
#include "acl/DataTypes/aclArray.h"
#include "acl/DataTypes/aclSubvector.h"
#include "acl/DataTypes/aclLocalArray.h"
#include "acl/Operators/aclElementFor.h"
#include "acl/Operators/aclElementIfElse.h"
#include "acl/Operators/aclElementExcerpt.h"
#include "acl/Kernels/aclKernel.h"
#include "aslUtilities.h"
#include <math.h>
#include <initializer_list>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>


using namespace acl;
using namespace std;

bool testCopy()
{
	cout << "Test of \"copy\" function..." << flush;
	Element vec0(new Array<cl_float> (10));
	
	vector<cl_float> input(10, 3);
	vector<cl_float> output(10, 1);

	copy(input, vec0);
	copy(vec0, output);

	bool status(output[3] == 3);
	errorMessage(status);

	return status;		
}


bool testKernel()
{
	cout << "Test of Kernel with double..." << flush;

	Element vec0(new Array<cl_double>(10));
	Element vec1(new Array<cl_double>(10));
	Element vec2(new Array<cl_double>(10));
	Element c(new Constant<cl_double>(2.));
	Element ind(new Index());


	Kernel k;
	{ 
		using namespace elementOperators;
		k.addExpression(operatorAssignment(vec2, c ));
		k.addExpression(operatorAssignment(vec0, c + powI(vec2, 3)));
		k.addExpression(operatorAssignment(vec1, ind));
	}
	k.setup();
	k.compute();

	vector<cl_double> output0(10), output1(10);
	copy(vec0, output0);	
	copy(vec1, output1);

	bool status(output0[9]<10.1 && output1[2]>2-1e-4 && output1[3]<3+1e-4);
	errorMessage(status);

	return status;		
}


bool testKernelSIMD()
{
	cout << "Test of KernelSIMD..." << flush;

	Element vec0(new Array<cl_float>(11));
	Element vec1(new Array<cl_float>(11));
	
	vector<cl_float> input0(11, 3);
	vector<cl_float> input1(11, 5);
	vector<cl_float> output(11, 0);
	vector<cl_float> expected({8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8});
	copy(input0, vec0);
	copy(input1, vec1);

	Kernel k(KERNEL_SIMD);
	{
		using namespace elementOperators;
		k.addExpression(operatorAssignmentSafe(vec1, vec0 + vec1));
	}

	k.setup();

	k.compute();
	copy(vec1, output);

	bool status(output == expected);
	errorMessage(status);

	return status;
}


bool testKernelSIMDUA()
{
	cout << "Test of KernelSIMDUA..." << flush;

	Element vec0(new Array<cl_float> (11));
	Element vec1(new Array<cl_float> (11));
	
	vector<cl_float> input0(11, 3);
	vector<cl_float> input1(11, 5);
	vector<cl_float> output(11, 0);
	vector<cl_float> expected({8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8});
	copy(input0, vec0);
	copy(input1, vec1);

	KernelConfiguration kConf(KERNEL_SIMDUA);
//	kConf.extensions.push_back("cl_amd_printf");
	Kernel k(kConf);
	{
		using namespace elementOperators;
		k.addExpression(operatorAssignmentSafe(vec1, vec0 + vec1));
//		k.addExpression(printfFunction("\"index: %d\\n\", index"));
	}

	k.setup();

	k.compute();
	copy(vec1, output);

	bool status(output == expected);
	errorMessage(status);

	return status;
}


bool testPrivateVariable()
{
	cout << "Test of kernel with PrivateVariable..." << flush;

	Element vec0(new Array<cl_float>(10));
	Element vec1(new Array<cl_float>(10));
	Element loc(new PrivateVariable<cl_float>());
	
	vector<cl_float> input1(10, 3);
	vector<cl_float> input2(10, 5);
	vector<cl_float> output(10, 1);

	copy(input1, vec0);
	copy(input2, vec1);

	Kernel k;
	{
		using namespace elementOperators;
		k.addExpression(operatorAssignment(loc, vec0 + vec1));
		k.addExpression(operatorAssignment(vec1, vec0 - vec1));
		k.addExpression(operatorAssignment(vec0, loc));
	}
	k.setup();
	
	k.compute();
	copy(vec0, output);


	bool status(output[2] ==8.);
	errorMessage(status);

	return status;
}


bool testPrivateArray()
{
	cout << "Test of kernel with PrivateArray..." << flush;
	
	vector<cl_int> inputGaIn({0, 4, 5});
	vector<cl_float> inputGaOut(3, 0);
	vector<cl_float> inputPa({-9, 2, 0, 15, 1, 3});
	vector<cl_float> output(3);
	vector<cl_float> expected({-9, 1, 3});

	Element gaIn(new Array<cl_int>(3));
	Element gaOut(new Array<cl_float>(3));
	Element pa(new PrivateArray<cl_float>(inputPa));
	shared_ptr<ElementExcerpt> ex(new ElementExcerpt(pa, gaIn));
	
	copy(inputGaIn, gaIn);


	Kernel k;
	{
		using namespace elementOperators;
		k.addExpression(operatorAssignment(gaOut, ex));
	}
	k.setup();

	k.compute();
	copy(gaOut, output);


	bool status(output == expected);
	errorMessage(status);

	return status;
}


bool testVariable()
{
	cout << "Test of Variable functionality..." << flush;

	Element vec0(new Array<cl_float> (10));
	shared_ptr<Variable<cl_float> > a(new Variable<cl_float> (1.));	
	
	vector<cl_float> output(10, 1);

	Kernel k;

	k.addExpression(elementOperators::operatorAssignment(vec0, a));	
	k.setup();

	k.compute();
	a->setValue(10.);
	k.compute();
	copy(vec0, output);
	
	bool status(output[2] ==10.);
	errorMessage(status);

	return status;		
}

bool testVariableReference()
{
	cout << "Test of VariableReference functionality..." << flush;

	Element vec0(new Array<cl_float> (10));
	float v(1.);
	Element a(new VariableReference<cl_float> (v));	
	
	vector<cl_float> output(10, 1);

	Kernel k;

	k.addExpression(elementOperators::operatorAssignment(vec0, a));	
	k.setup();

	k.compute();
	v=10.;
	k.compute();
	copy(vec0, output);

	bool status(output[2] ==10.);
	errorMessage(status);

	return status;		
}


bool testSelect()
{
	cout << "Test of select function..." << flush;

	Element vec0(new Array<cl_double> (10));
	Element c0(new Constant<cl_double> (2.1));

	vector<cl_double> input(10, 3.);
	vector<cl_double> output(10, 1.);

	copy(input,vec0);

	Kernel k;
	{
		using namespace elementOperators;
		k.addExpression(operatorAssignment(vec0, 
		                                   select(vec0-c0, 
		                                          vec0*vec0, 
		                                          convert(TYPE_SELECT[TYPE_DOUBLE],
		                                                  vec0 > c0,
		                                                  false))));
	}	
	k.setup();

	k.compute();
	copy(vec0, output);
	
	bool status(output[2] ==9.);
	errorMessage(status);

	return status;		
}


bool testSubvector()
{
	cout << "Test of Subvector..." << flush;
	cl_float init[] = {16, 2, 77, 29, 23, 16, 2, 77, 29, 23};
	shared_ptr<Array<cl_float> > vec0(new Array<cl_float>(10));
	
	vector<cl_float> input(init, init + sizeof(init) / sizeof(cl_float) );
	vector<cl_float> output(2);

	copy(input, vec0);
	Element subvec0(new Subvector<cl_float> (vec0, 2, 0));
	copy(subvec0, output);

	bool status(output[0]==16);
	errorMessage(status);

	return status;
}


bool testSwapBuffers()
{
	cout << "Test of Swap functionality..." << flush;
	shared_ptr<Array<cl_float> > vec0(new Array<cl_float>(10));
	shared_ptr<Array<cl_float> > vec1(new Array<cl_float>(10));
	
	vector<cl_float> input0(10, 1);
	vector<cl_float> input1(10, 2);
	vector<cl_float> output(10, 10);

	copy(input0, vec0);
	copy(input1, vec1);
	swapBuffers(*vec1,*vec0);
	copy(vec0, output);

	bool status(output[3] == 2);
	errorMessage(status);

	return status;	
}


bool testLocalArray()
{
	cout << "Test of LocalArray and syncCopy with barrier()..." << flush;

	KernelConfiguration kConf(KERNEL_BASIC);
	kConf.local = true;

	unsigned int groupsNumber = 5;
	unsigned int groupSize = 2;
	
	Element vec0(new Array<cl_float>(groupSize * groupsNumber));
	Element vec1(new Array<cl_float>(groupSize * groupsNumber));
	Element loc0(new LocalArray<cl_float>(groupSize));
	Element loc1(new LocalArray<cl_float>(groupSize));
	Element groupID(new GroupID());
	Element cGroupSize(new Constant<cl_uint>(groupSize));
	Element c0(new Constant<cl_uint>(0));
	
	vector<cl_float> input0(groupSize * groupsNumber, 3);
	vector<cl_float> input1(groupSize * groupsNumber, 5);
	vector<cl_float> output(groupSize * groupsNumber, 0);
	vector<cl_float> expected({2, 2, 2, 2, 2, 2, 2, 2, 2, 2});
	copy(input0, vec0);
	copy(input1, vec1);

	Kernel k(kConf);
	k.setGroupsNumber(groupsNumber);
	{
		using namespace elementOperators;
		k.addExpression(syncCopy(vec0, loc0, cGroupSize * groupID, c0, cGroupSize));
		k.addExpression(syncCopy(vec1, loc1, cGroupSize * groupID, c0, cGroupSize));
		k.addExpression(barrier());
		k.addExpression(operatorAssignment(loc1, loc1 - loc0));
		k.addExpression(barrier("CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE"));
		k.addExpression(syncCopy(loc1, vec1, c0, cGroupSize * groupID, cGroupSize));
	}
	k.setup();

	
	k.compute();
	copy(vec1, output);

	bool status(output == expected);
	errorMessage(status);

	return status;

}


int main()
{
	bool allTestsPassed(true);

	allTestsPassed &= testCopy();
	allTestsPassed &= testKernel();
	allTestsPassed &= testKernelSIMD();
	allTestsPassed &= testKernelSIMDUA();
	allTestsPassed &= testPrivateVariable();
	allTestsPassed &= testPrivateArray();
	allTestsPassed &= testVariable();
	allTestsPassed &= testVariableReference();
	allTestsPassed &= testSelect();
	allTestsPassed &= testSwapBuffers();
	allTestsPassed &= testLocalArray();
	allTestsPassed &= testSubvector();
	
	return allTestsPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}