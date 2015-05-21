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


#include "aclOperatorBinary.h"
#include "../aclUtilities.h"
#include <algorithm>

using namespace std;
using namespace asl;

namespace acl
{

	OperatorBinary::OperatorBinary(Element a1, Element a2, const string & description):
		// first operand dictates the typeID of the binary operator
		ElementBase(false, 0, a1->getTypeID())
	{
		if ( compatible(a1, a2) )
		{
			size = max(a1->getSize(), a2->getSize());

			if (a1->getQueue().get() != 0)
				queue = a1->getQueue();
			else
				queue = a2->getQueue();
			
			e1 = a1;
			e2 = a2;
		}
		else
		{
			errorMessage("Operands of binary operator (" + description + ") are incompatible. \
				 		 Either they reside on different devices or their sizes do not match: "
			             + numToStr(a1->getSize()) + " and " + numToStr(a2->getSize()));
		}
	}
	
	void OperatorBinary::addToKernelSource(vector<Element> & arguments,
			                               vector<Element> & localDeclarations) const
	{
		addElementToKernelSource(e1, arguments, localDeclarations);
		addElementToKernelSource(e2, arguments, localDeclarations);
	}


	void OperatorBinary::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
	{
	}


	string OperatorBinary::getName() const
	{
		return "";
	}


	string OperatorBinary::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	string OperatorBinary::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}

} // namespace acl
