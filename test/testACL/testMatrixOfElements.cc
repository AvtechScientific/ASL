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
	\example testMatrixOfElements.cc
 */

#include "acl/acl.h"
#include "acl/aclGenerators.h"
#include "acl/aclMath/aclMatrixOfElements.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include "math/aslMatrices.h"
#include "acl/Kernels/aclKernel.h"
#include "acl/DataTypes/aclArray.h"

using namespace asl;
using namespace acl;
using namespace std;

bool testMatrixOperations()
{
	cout << "Test of \"Matrix Operations\" function..." << flush;

	VectorOfElements vec0(3);
	VectorOfElements vec1(1);
	copy(generateVEData<cl_float>(10u,3u),vec0);
	copy(generateVEData<cl_float>(10u,1u),vec1);

	Kernel k;
	{ 
		using namespace elementOperators;
		k << (vec0=generateVEConstant(0.1f,1.f,2.f));
		k << (vec1=(elementProduct(generateVEConstant(0.f,1.f,2.f),vec0)*generateVEConstant(1.f,0.f,2.f))*vec0);
	}
	k.setup();
	k.compute();

	vector<cl_float> output(10);
	copy(vec1[0], output);
	
	bool status(output[1] == 20.5);
	errorMessage(status);

	return status;
}


bool testSystemSolve()
{
	cout << "Test of \"System Solve Cramer's rule\" function..." << flush;

	VectorOfElements vecB(2);
	VectorOfElements vecX(2);
	copy(generateVEData<cl_float>(10u,2u),vecB);
	copy(generateVEData<cl_float>(10u,2u),vecX);
	
	auto matA(generateMEConstant(makeAMatr(makeAVec(4.,1.),makeAVec(1.,3.))));

	Kernel k;
	{ 
		using namespace elementOperators;
		k << (vecB=generateVEConstant(1.f,2.f));
		k << gcSolveSystem(matA,vecB,vecX);
	}
	k.setup();
	k.compute();

	vector<cl_float> output(10);
	copy(vecX[0], output);
	
	bool status(output[1] > 0.09 && output[1] < .1);
	errorMessage(status);

	return status;
}

bool testSystemSolveCG()
{
	cout << "Test of \"System Solve congugate gradient method\" function..." << flush;

	VectorOfElements vecB(2);
	VectorOfElements vecX(2);
	copy(generateVEData<cl_float>(10u,2u),vecB);
	copy(generateVEData<cl_float>(10u,2u),vecX);
	
	auto matA(generateMEConstant(makeAMatr(makeAVec(4.,1.),makeAVec(1.,3.))));

	Kernel k;
	{ 
		using namespace elementOperators;
		k << (vecB=generateVEConstant(1.f,2.f));
		k << gcSolveSystemCG(matA,vecB,vecX);
	}
	k.setup();
	k.compute();

	vector<cl_float> output(10);
	copy(vecX[0], output);
	
	bool status(output[1] > 0.09 && output[1] < .1);
	errorMessage(status);

	return status;		
}


int main()
{
	bool allTestsPassed(true);

	allTestsPassed &= testMatrixOperations();
	allTestsPassed &= testSystemSolve();
	allTestsPassed &= testSystemSolveCG();

	return allTestsPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}