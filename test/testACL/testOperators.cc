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
	\example testOperators.cc
 */


#include "acl/acl.h"
#include "acl/DataTypes/aclIndex.h"
#include "acl/DataTypes/aclConstant.h"
#include "acl/DataTypes/aclVariable.h"
#include "acl/DataTypes/aclVariableReference.h"
#include "acl/DataTypes/aclPrivateVariable.h"
#include "acl/DataTypes/aclArray.h"
#include "acl/DataTypes/aclLocalArray.h"
#include "acl/DataTypes/aclSubvector.h"
#include "acl/Operators/aclElementFor.h"
#include "acl/Operators/aclElementIfElse.h"
#include "acl/Operators/aclElementParser.h"
#include "acl/Kernels/aclKernel.h"
#include "aslUtilities.h"
#include <math.h>
#include <initializer_list>

using namespace acl;
using namespace std;


bool testIfElse()
{
	cout << "Test of If-Else..." << flush;

	using namespace elementOperators;
	shared_ptr<Variable<cl_int> > a(new Variable<cl_int> (15));
	Element c0(new Constant<cl_int> (15));
	Element c1(new Constant<cl_int> (3));
	Element vec(new Array<cl_float> (11));
	vector<cl_float> input(11, 2); 
	vector<cl_float> output(11, 0);
	vector<cl_float> expected({8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8});
	copy(input, vec);

	// Test if
	shared_ptr<ElementIfElse> ifElse_TestIf(new ElementIfElse(isEqual(a, c0)));
	ifElse_TestIf->addBodyExpressionIf(operatorAssignment(vec, vec + c1));
	ifElse_TestIf->addBodyExpressionElse(operatorAssignment(vec, vec - c1));

	// Test else
	shared_ptr<ElementIfElse> ifElse_TestElse(new ElementIfElse(isEqual(a, c1)));
	ifElse_TestElse->addBodyExpressionIf(operatorAssignment(vec, vec - c1));
	ifElse_TestElse->addBodyExpressionElse(operatorAssignment(vec, vec + c1));


	Kernel k;

	k.addExpression(ifElse_TestIf);
	k.addExpression(ifElse_TestElse);
	k.setup();
	k.compute();
	copy(vec, output);

	bool status(output == expected);
	errorMessage(status);

	return status;		
}


bool testParser()
{
	cout << "Test of Parser..." << flush;

	using namespace elementOperators;
	shared_ptr<Variable<cl_int> > a(new Variable<cl_int> (15));
	Element c0(new Constant<cl_int> (15));
	Element c1(new Constant<cl_int> (3));
	Element vec(new Array<cl_float> (11));
	vector<cl_float> input(11, 2);
	vector<cl_float> output(11, 0);
	vector<cl_float> expected({35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35});

	copy(input, vec);

	string statement("a + c0 + c1 + vec");
	shared_ptr<ElementParser> parser(new ElementParser());
	parser->addElementNamePair(a, "a");
	parser->addElementNamePair(c0, "c0");
	parser->addElementNamePair(c1, "c1");
	parser->addElementNamePair(vec, "vec");
	parser->setStatement(statement);
	
	Kernel k;

	k.addExpression(operatorAssignment(vec, parser));
	k.setup();
	k.compute();

	copy(vec, output);

	bool status(output == expected);
	errorMessage(status);

	return status;		
}

bool testAtomicSum()
{
	cout << "Test of Atomic Sum..." << flush;

	using namespace elementOperators;
	Element c(new Constant<cl_int> (6));
	Element vec(new Array<cl_int> (11));
	vector<cl_int> input(11, 2); 
	vector<cl_int> output(11, 0);
	vector<cl_int> expected({8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8});
	copy(input, vec);

	KernelConfiguration kConf(KERNEL_BASIC);
	kConf.extensions.push_back("cl_khr_global_int32_base_atomics");
	Kernel k(kConf);

	k.addExpression(atomic_add(vec, c));
	k.setup();
	k.compute();
	copy(vec, output);

	bool status(output == expected);
	errorMessage(status);

	return status;		
}

int main()
{
	bool allTestsPassed(true);

	allTestsPassed &= testIfElse();
	allTestsPassed &= testParser();
	allTestsPassed &= testAtomicSum();

	return allTestsPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}
