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


#include "aslSingleKernelNM.h"
#include <data/aslDataWithGhostNodes.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <acl/Kernels/aclKernel.h>
#include <acl/acl.h>


namespace asl
{
	SingleKernelNM::SingleKernelNM(const acl::KernelConfiguration & kernelConfig):
		kernel(new acl::Kernel(kernelConfig))
	{
	}
	
	SingleKernelNM::~SingleKernelNM()
	{
	}

	void SingleKernelNM::init()
	{
		init0();
		kernel->setup();
	}

	void SingleKernelNM::execute()
	{
		preProcessing();
		kernel->compute();
		postProcessing();
	}

	void SingleKernelNM::preProcessing()
	{
	}

	void SingleKernelNM::postProcessing()
	{
	}

	SingleKernelMapNM::SingleKernelMapNM(const acl::KernelConfiguration & kernelCongig):
		SingleKernelNM(kernelCongig)
	{
	}

	SingleKernelMapNM::SingleKernelMapNM(Field map, const acl::KernelConfiguration & kernelCongig):
		SingleKernelNM(kernelCongig),
		map(map)			
	{
	}

	void SingleKernelMapNM::setMap(Field m)
	{
		map=m;
	}

	void SingleKernelMapNM::initMapInfrastructure(acl::ExpressionContainer & k)
	{
		if (map.get()!=0)
		{
			k.addExpression(acl::elementOperators::
				                ifElse(acl::elementOperators::any((map->getEContainer() <= 0.)[0]),
				                       {acl::elementOperators::returnStatement()}, 
				                       {} ));
		}
	}

	SingleKernelMapNM::~SingleKernelMapNM()
	{
	}

} // asl

