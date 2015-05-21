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


#include "aslDataUtilities.h"
#include "aslGenerators.h"
#include "acl/aclElementBase.h"
#include "math/aslVectorsDynamicLengthOperations.h"
#include "acl/aclGenerators.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include "acl/DataTypes/aclIndex.h"
#include "acl/DataTypes/aclConstant.h"
#include "acl/acl.h"
#include "acl/Kernels/aclKernel.h"
#include "math/aslIndex2Position.h"
#include "acl/aclUtilities.h"
#include "acl/Operators/aclElementIfElse.h"


namespace asl
{

	void checkDimensionsNumCompatibility(AbstractData & source,
	                                     const AVec<int> & size)
	{
		if (nD(source.getBlock()) !=  nD(size))
			errorMessage("uploadToLocalMem() - dimensions mismatch");
	}


	void checkTessellability(AbstractData & source,
	                         const AVec<int> & size)
	{
		const AVec<int> s(source.getBlock().getSize());
		for (unsigned int i(0); i < nD(s); ++i)
			if (s[i] % size[i] != 0)
				errorMessage("uploadToLocalMem() - size of source is not evenly divisible by the block\'s size");
	}
	
	
	acl::VectorOfElements uploadToLocalMem(AbstractData & source,
	                                       const AVec<int> & size,
	                                       unsigned int groupSize,
	                                       acl::Kernel & kernel)
	{
		checkDimensionsNumCompatibility(source,size);
		checkTessellability(source, size);
		
		AVec<int> groupTableSize(divisionOfElements(source.getBlock().getSize(),
		                                            size));

		unsigned int totalPointsNum = productOfElements(size);
		unsigned int uploadsPerWorkItem = totalPointsNum / groupSize;
		unsigned int uploadsPerWorkItemRemainder = totalPointsNum % groupSize;

		unsigned int componentsNum = source.getEContainer().size();
		acl::VectorOfElements uploadedData(componentsNum);

		copy(generateVELocalArray(totalPointsNum,
		                          source.getEContainer()[0]->getTypeID(),
		                          componentsNum),
		     uploadedData);

		Block localBlock(size);
		Index2PositionDiscreteACL ind(localBlock, false);
		Block groupTableBlock(groupTableSize);
		Index2PositionDiscreteACL indGroup(groupTableBlock);

		acl::VectorOfElements transformVec(acl::generateVEConstant(source.getBlock().c2iTransformVector));	

		acl::VectorOfElements offset(productOfElements(acl::excerpt(indGroup.positionWithInit,
		                                                            acl::generateVEGroupID()),
		                                               acl::generateVEConstant(size)));

		acl::VectorOfElements singleUpload(uploadedData = acl::excerpt(source.getEContainer(),
		                                                               (ind.positionWithInit + offset) * transformVec));

		acl::VectorOfElements index(acl::generateVEIndex(groupSize));

		for (unsigned int i = 0; i < uploadsPerWorkItem; ++i)
		{
			kernel << (excerpt(singleUpload,
			                   groupSize * i + index));
		}

		acl::Element uploadsPerWorkItemRem(new acl::Constant<cl_int>(uploadsPerWorkItemRemainder));
		acl::Element indexRemainder(new acl::Index());

		using namespace acl::elementOperators;
		acl::Element ife(ifElse(indexRemainder < uploadsPerWorkItemRem,
		                        excerpt(singleUpload,
		                                groupSize * uploadsPerWorkItem + index),
		                        {}));
		kernel.addExpression(ife);

		return uploadedData;
	}

	SPDataWrapperACL generateSubData(SPDataWrapperACL d, AVec<int> a, AVec<int> b)
	{
		const Block & bIn(d->getBlock());
		unsigned int nd(nD(a));
		Block bOut(b - a + AVec<int>(nd,1), bIn.dx, bIn.position+bIn.dx*AVec<>(a));

		acl::VectorOfElements c(d->getEContainer().size());

		Index2PositionDiscreteACL t2p(bOut);
		acl::VectorOfElements newIndex((t2p.positionWithInit + a) *
		                                bIn.c2iTransformVector);
		
		copy(acl::excerpt(d->getEContainer(), newIndex), c);		

		return generateDataContainer_SP(bOut,c);
	}


} // namespace asl
