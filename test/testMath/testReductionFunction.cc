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
#include <acl/aclUtilities.h>
#include "acl/aclMath/aclReductionAlgGenerator.h"
#include "acl/aclGenerators.h"
#include "aslUtilities.h"
#include "acl/aclMath/aclVectorOfElements.h"

using namespace acl;

void testSum()
{
	unsigned int n(101);
	auto v(generateVEData<float>(n,1u));
	initData(v, generateVEConstant(2));
	auto summator(generateSumAlg<float>(v));
	summator->generateAlg();
	summator->compute();
	bool b(asl::approxEqual(summator->res.v()[0],2.f*n));
	cout<<"testSum: "<<(b? "Ok": "Error ")<<endl;
}

void testSum1()
{
	unsigned int n(100001);
	VectorOfElements v1(generateVEData<float>(n,1u));
	VectorOfElements v2(generateVEData<float>(n,1u));
	initData(v1, generateVEConstant(2));
	initData(v2, generateVEConstant(3));
	auto summator(generateSumAlg<float>(v1*v2));
	summator->generateAlg();
	summator->compute();
	bool b(asl::approxEqual(summator->res.v()[0],6.f*n));
	cout<<"testSum1: "<<(b? "Ok": "Error")<<endl;
}

bool testMin()
{
	VectorOfElements vI(generateVEIndex());
	VectorOfElements v1(generateVEData<float>(101u,1u));
	initData(v1, generateVEConstant(2));
	auto minimizer(generateMinAlg<float>(v1*((vI-100)*(vI-100)+3)));
	minimizer->generateAlg();
	minimizer->compute();
	bool b(asl::approxEqual(minimizer->res.v()[0],6.f));
	cout<<"testMin: "<<(b? "Ok": "Error")<<endl;
	return b;
}

bool testMax()
{
	VectorOfElements vI(generateVEIndex());
	VectorOfElements v1(generateVEData<float>(100001u,1u));
	initData(v1, generateVEConstant(2));
	auto maximizer(generateMaxAlg<float>(v1*((1000.-vI)*(vI-1000.)-10.)));
	maximizer->generateAlg();
	maximizer->compute();
	bool b(asl::approxEqual(maximizer->res.v()[0],-20.f));
	cout<<"testMax: "<<(b? "Ok": "Error")<<endl;
	return b;
}

bool testProduct()
{
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
	bool b(asl::approxEqual(alg->res.v()[0],256));
	cout<<"testProduct: "<<(b? "Ok": "Error")<<endl;
	return b;
}


int main()
{
	testSum();
	testSum1();
	testMin();
	testMax();
	testProduct();
	return 0;
}
