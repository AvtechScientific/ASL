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


#include "aclKernelMerger.h"
#include "aclKernel.h"
#include <acl/acl.h>
#include "acl/Operators/aclElementIfElse.h"
#include "acl/DataTypes/aclIndex.h"
#include "acl/DataTypes/aclConstant.h"
#include <aslUtilities.h>

namespace acl
{
	KernelMerger::KernelMerger():
		kernel(new Kernel())
	{
	}


	void KernelMerger::compute()
	{
		kernel->compute();		
	}

	unsigned int KernelMerger::getKernelSize(unsigned int i)
	{
		const KernelConfiguration &kc(kernels[i]->kernelConfig);
		unsigned int size(kernels[i]->getSize());
		if (kc.vectorWidth > 1)
		{
			size += (kc.vectorWidth - size % kc.vectorWidth) % kc.vectorWidth;
				return size;
		}
		return size;
	}
		

	std::string KernelMerger::getKernelSource()
	{
		return kernel->getKernelSource();
	}
	

	void KernelMerger::clear()
	{
		kernels.clear();
		offsets.clear();
		kernel->clear();
	}


	void KernelMerger::addKernel(SPKernel k)
	{
		kernels.push_back(k);
	}


	void KernelMerger::addKernel(const KernelMerger & km)
	{
		for(unsigned int i(0); i < km.kernels.size(); ++i)
		kernels.push_back(km.kernels[i]);
	}
		

	unsigned int KernelMerger::getSize()
	{
		return size;
	}


	void KernelMerger::computeOffsets()
	{
		size = getKernelSize(0);
		offsets.resize(kernels.size() - 1);
		for (unsigned int i(1); i < kernels.size(); ++i)
		{
			const KernelConfiguration &kc(kernels[0]->kernelConfig);
			offsets[i-1] = kc.unaligned ? size : size / kc.vectorWidth;
			size += getKernelSize(i);
		}
	}


	void KernelMerger::setup()
	{
		if (!kernels.size())
			asl::errorMessage("KernelMerger::setup() : no kernels have been added."); 
		if (kernels.size() > 1)
		{
			kernel.reset(new Kernel(kernels[0]->kernelConfig));
			computeOffsets();
			kernel->addExpression(castSpliter(0, kernels.size() - 1));		
		}
		else
			kernel = kernels[0];
		kernel->setup();
	}


	void addToIfBody(std::shared_ptr<ElementIfElse> ife,
	                 std::shared_ptr<Kernel> k,
	                 unsigned int size)
	{
		Element ind(new Index(size));
		for (unsigned int i(0); i < k->expression.size(); ++i)
			ife->addBodyExpressionIf(elementOperators::excerpt(k->expression[i], ind));
	}
	

	void addToElseBody(std::shared_ptr<ElementIfElse> ife,
	                   std::shared_ptr<Kernel> k,
	                   unsigned int size)
	{
		Element ind(new Index(size));
		for (unsigned int i(0); i < k->expression.size(); ++i)
			ife->addBodyExpressionElse(elementOperators::excerpt(k->expression[i],ind));
		
	}


	Element KernelMerger::castSpliter(unsigned int i1, unsigned int i2)
	{
		Element ind(new Index(size));
		if (i2 - i1 == 1)
		{
			using namespace elementOperators;
			Element c(new Constant<int>(offsets[i1]));
			shared_ptr<ElementIfElse> ife(new ElementIfElse(ind < c));
			if (i1 > 0)
			{
				Element c(new Constant<int>(offsets[i1 - 1]));
				ife->addBodyExpressionIf(ind -= c);
			}
			addToIfBody(ife, kernels[i1], size);
			ife->addBodyExpressionElse(ind -= c);
			addToElseBody(ife, kernels[i2], size);
			return ife;
		}
		if (i2 - i1 == 2)
		{
			using namespace elementOperators;
			Element c(new Constant<int>(offsets[i1]));
			shared_ptr<ElementIfElse> ife(new ElementIfElse(ind < c));
			if (i1 > 0)
			{
				Element c(new Constant<int>(offsets[i1 - 1]));
				ife->addBodyExpressionIf(ind -= c);
			}
			addToIfBody(ife, kernels[i1],size);
			ife->addBodyExpressionElse(castSpliter(i1 + 1, i2));
			return ife;
		}

		using namespace elementOperators;
		unsigned int midi(i1 + (i2 - i1) / 2);
		Element c(new Constant<int>(offsets[midi]));
		shared_ptr<ElementIfElse> ife(new ElementIfElse(ind < c));
		ife->addBodyExpressionIf(castSpliter(i1, midi));
		ife->addBodyExpressionElse(castSpliter(midi + 1, i2));

		return ife;
	}	
	
} // namespace acl
