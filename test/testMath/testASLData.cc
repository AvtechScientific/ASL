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
	\example testASLData.cc
 */

#include "acl/acl.h"
#include "acl/aclGenerators.h"
#include "acl/DataTypes/aclGroupID.h"
#include "acl/DataTypes/aclConstant.h"
#include "aslGenerators.h"
#include "aslUtilities.h"
#include "acl/Kernels/aclKernel.h"
#include "data/aslDataUtilities.h"
#include <acl/aclMath/aclVectorOfElements.h>
#include <data/aslDataWithGhostNodes.h>

using namespace acl;
using namespace asl;
using namespace std;

bool testSimpleKernel()
{
	Block bl(makeAVec(10,5),0.1);
	auto a(generateDataContainerACL_SP<double>(bl,1,1));
	auto b(generateDataContainerACL_SP<double>(bl,1,1));

	Kernel k;
	
	k<<(a->getEContainer()=acl::generateVEConstant<double>(0));
	k<<(b->getEContainer()=acl::generateVEConstant<double>(0));

	k.setup();
	cout<<k.getKernelSource ()<<endl;

	return true;
}


bool testInitData()
{
	Block bl(makeAVec(10,5),0.1);
	auto a(generateDataContainerACL_SP<double>(bl, 1, 1));

	initData(a->getEContainer(),acl::generateVEConstant<double>(0));

	return true;
}


bool testUploadToLocalMem()
{
	cout << "Test of UploadToLocalMem()..." << flush;

	unsigned int componentsNum = 2;
	unsigned int groupSize = 27;
	KernelConfiguration kConf(KERNEL_BASIC);
	kConf.local = true;
	Kernel kernel(kConf);
	Element groupID(new GroupID());
	Element c0(new Constant<cl_uint>(0));

	asl::AVec<int> totalDimensions(asl::makeAVec(10, 10, 10));
	// dx = 1
	asl::Block block(totalDimensions, 1);

	// boundary = 0
	auto source(asl::generateDataContainerACL_SP<float>(block,
	                                                     componentsNum,
	                                                     0u));
	// initialize source with value "13"
	acl::initData(source->getEContainer(), acl::generateVEConstantN(componentsNum,
	                                                                13));

	// boundary = 0
	auto destination(asl::generateDataContainerACL_SP<float>(block,
	                                                          componentsNum,
	                                                          0u));
	// initialize source with value "27"
	acl::initData(destination->getEContainer(), acl::generateVEConstantN(componentsNum,
	                                                                     27));


	asl::AVec<int> portionDimensions(asl::makeAVec(5, 5, 5));

	unsigned int portionSize = productOfElements(portionDimensions);
	unsigned int totalSize = productOfElements(totalDimensions);
	unsigned int groupsNumber = totalSize / portionSize;
	Element cPortionSize(new Constant<cl_uint>(portionSize));
	kernel.setGroupsNumber(groupsNumber);

	VectorOfElements localSource(componentsNum);
	copy(uploadToLocalMem(*source, portionDimensions, groupSize, kernel), localSource);

	using namespace elementOperators;
	kernel.addExpression(barrier("CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE"));
	// overwrite destination with localSource
	for (unsigned int i = 0; i < componentsNum; ++i)
		kernel.addExpression(syncCopy(localSource[i], destination->getEContainer()[i],
		                              c0, cPortionSize * groupID, cPortionSize));

	kernel.setup();
	kernel.compute();


	bool status(true);
	vector<cl_float> src(totalSize);
	vector<cl_float> dst(totalSize);
	for (unsigned int i = 0; i < componentsNum; ++i)
	{
		copy(source->getEContainer()[i], src);
		copy(destination->getEContainer()[i], dst);

		cout << src << endl << endl << dst << endl;
		
		status &= (src == dst);
	}
	
	errorMessage(status);

	return status;
}

int main()
{
	bool allTestsPassed(true);

	allTestsPassed &= testSimpleKernel();
	allTestsPassed &= testInitData();
	allTestsPassed &= testUploadToLocalMem();
		
	return allTestsPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}