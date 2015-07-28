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
	\example testReductionFunction.cc
 */

#include "acl/Kernels/aclKernel.h"
#include "acl/aclUtilities.h"
#include "acl/aclMath/aclReductionAlgGenerator.h"
#include "acl/aclGenerators.h"
#include "aslUtilities.h"
#include "acl/aclMath/aclVectorOfElements.h"

using namespace acl;

bool testSum()
{
	cout << "testSum..." << flush;
	unsigned int n(101);
	auto v(generateVEData<float>(n,1u));
	initData(v, generateVEConstant(2));
	auto summator(generateSumAlg<float>(v));
	summator->generateAlg();
	summator->compute();
	bool status(asl::approxEqual(summator->res.v()[0],2.f*n));
	asl::errorMessage(status);

	return status;
}


bool testSum1()
{
	cout << "testSum1..." << flush;
	unsigned int n(100001);
	VectorOfElements v1(generateVEData<float>(n,1u));
	VectorOfElements v2(generateVEData<float>(n,1u));
	initData(v1, generateVEConstant(2));
	initData(v2, generateVEConstant(3));
	auto summator(generateSumAlg<float>(v1*v2));
	summator->generateAlg();
	summator->compute();
	bool status(asl::approxEqual(summator->res.v()[0],6.f*n));
	asl::errorMessage(status);

	return status;
}


bool testMin()
{
	cout << "testMin..." << flush;
	VectorOfElements vI(generateVEIndex());
	VectorOfElements v1(generateVEData<float>(101u,1u));
	initData(v1, generateVEConstant(2));
	auto minimizer(generateMinAlg<float>(v1*((vI-100)*(vI-100)+3)));
	minimizer->generateAlg();
	minimizer->compute();
	bool status(asl::approxEqual(minimizer->res.v()[0],6.f));
	asl::errorMessage(status);

	return status;
}


bool testMax()
{
	cout << "testMax..." << flush;
	VectorOfElements vI(generateVEIndex());
	VectorOfElements v1(generateVEData<float>(100001u,1u));
	initData(v1, generateVEConstant(2));
	auto maximizer(generateMaxAlg<float>(v1*((1000.-vI)*(vI-1000.)-10.)));
	maximizer->generateAlg();
	maximizer->compute();
	bool status(asl::approxEqual(maximizer->res.v()[0],-20.f));
	asl::errorMessage(status);

	return status;
}

bool testProduct()
{
	cout << "testProduct..." << flush;
	typedef double FT;
	VectorOfElements vI(generateVEIndex());
	VectorOfElements v1(generateVEData<FT>(100001u,1u));
	initData(v1, generateVEConstant(2));
	auto alg(generateProductAlg<FT>(select(generateVEConstant(1.),
	                                       v1,
	                                       vI >=1000 && vI <= 1007, 
	                                       acl::typeToTypeID<FT>())));
	alg->generateAlg();
	alg->compute();
	bool status(asl::approxEqual(alg->res.v()[0],256));
	asl::errorMessage(status);

	return status;
}


int main()
{
	bool allTestsPassed(true);

	allTestsPassed &= testSum();
	allTestsPassed &= testSum1();
	allTestsPassed &= testMin();
	allTestsPassed &= testMax();
	allTestsPassed &= testProduct();

	return allTestsPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}