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
	\example testKernelMerger.cc
 */


#include "acl/acl.h"
#include "acl/DataTypes/aclConstant.h"
#include "acl/DataTypes/aclArray.h"
#include "acl/Kernels/aclKernel.h"
#include "acl/Kernels/aclKernelMerger.h"
#include "aslUtilities.h"
#include <math.h>
#include <initializer_list>

using namespace acl;
using namespace std;

bool testKernelMerger()
{
	cout << "Test of \"KernelMerger\" functionality..." << flush;
	ElementData vec0(new Array<cl_float> (10));
	ElementData vec1(new Array<cl_float> (5));
	ElementData vec2(new Array<cl_float> (8));
	ElementData vec3(new Array<cl_float> (20));	

	Element c0(new Constant<cl_double>(2));
	Element c1(new Constant<cl_double>(1));
	Element c2(new Constant<cl_double>(4));
	Element c3(new Constant<cl_double>(7));
	
	SPKernel k0(new Kernel());
	SPKernel k1(new Kernel());
	SPKernel k2(new Kernel());
	SPKernel k3(new Kernel());
	{ 
		using namespace elementOperators;
		k0->addExpression(operatorAssignment (vec0, c0));
		k1->addExpression(operatorAssignment (vec1, c1));
		k2->addExpression(operatorAssignment (vec2, c2));
		k3->addExpression(operatorAssignment (vec3, c3));
	}

	KernelMerger km;
	km.addKernel(k0);
	km.addKernel(k1);
	km.addKernel(k2);
//	km.addKernel(k3);

	km.setup();
	km.compute();

	bool status((acl::map<float>(vec0).get()[9] == 2) && 
	            (acl::map<float>(vec1).get()[3] == 1) &&
	            (acl::map<float>(vec2).get()[7] == 4));// &&
//	            (acl::map<float>(vec3).get()[19] == 7));
	errorMessage(status);
	cout << km.getKernelSource() << endl;
	return status;		
}

int main()
{
	bool allTestsPassed(true);

	allTestsPassed &= testKernelMerger();
	
	return allTestsPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}