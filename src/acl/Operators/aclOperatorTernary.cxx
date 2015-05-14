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


#include "aclOperatorTernary.h"
#include "../aclUtilities.h"
#include <algorithm>

using namespace std;
using namespace asl;

namespace acl
{

	OperatorTernary::OperatorTernary(Element a1, Element a2, Element a3):
		// first operand dictates the typeID of the ternary operator
		ElementBase(false, 0, a1->getTypeID())
	{
		// if three elements have the same size or size of at least one of them is 0
		if (compatible(a1,a2) && compatible(a2,a3) && compatible(a1,a3)) 
		{
			size = max(max(a1->getSize(), a2->getSize()), a3->getSize());

			if (a1->getQueue().get() != 0)
				queue = a1->getQueue();
			else
				if (a2->getQueue().get() != 0)
					queue = a2->getQueue();
				else
					queue = a3->getQueue();
			
			e1 = a1;
			e2 = a2;
			e3 = a3;
		}
		else
		{
			errorMessage("Operands of ternary operator are incompatible. \
				 		 Either they reside on different devices or their sizes do not match: "
			             + numToStr(a1->getSize()) + " and "
			             + numToStr(a2->getSize()) + " and "
			             + numToStr(a3->getSize()));
		}
	}


	void OperatorTernary::addToKernelSource(vector<Element> & arguments,
	                                        vector<Element> & localDeclarations) const
	{
		addElementToKernelSource(e1, arguments, localDeclarations);
		addElementToKernelSource(e2, arguments, localDeclarations);
		addElementToKernelSource(e3, arguments, localDeclarations);
	}


	void OperatorTernary::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
	{
	}


	string OperatorTernary::getName() const
	{
		return "";
	}


	string OperatorTernary::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	string OperatorTernary::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}

} // namespace acl
