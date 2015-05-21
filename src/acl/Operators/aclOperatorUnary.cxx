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


#include "aclOperatorUnary.h"
#include "../aclUtilities.h"

using namespace std;

namespace acl
{

	OperatorUnary::OperatorUnary(Element a):
		ElementBase(false, 0, a->getTypeID()),
		e(a)
	{
		size = a->getSize();
		queue = a->getQueue();
	}

	OperatorUnary::OperatorUnary(Element a, TypeID type):
		ElementBase(false, 0, type),
		e(a)
	{
		size = a->getSize();
		queue = a->getQueue();
	}

	void OperatorUnary::addToKernelSource(vector<Element> & arguments,
			                              vector<Element> & localDeclarations) const
	{
		addElementToKernelSource(e, arguments, localDeclarations);
	}


	void OperatorUnary::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
	{
	}


	string OperatorUnary::getName() const
	{
		return "";
	}


	string OperatorUnary::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	string OperatorUnary::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


} // namespace acl
