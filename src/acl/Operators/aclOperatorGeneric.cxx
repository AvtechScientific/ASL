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


#include "aclOperatorGeneric.h"

using namespace std;

namespace acl
{

	OperatorGeneric::OperatorGeneric(const string & operatorStr_):
		// typeID fictitious
		ElementBase(false, 0, TYPE_INT),
		operatorStr(operatorStr_)
	{
	}


	string OperatorGeneric::str(const KernelConfiguration & kernelConfig) const
	{
		return operatorStr;
	}


	void OperatorGeneric::addToKernelSource(vector<Element> & arguments,
	                                        vector<Element> & localDeclarations) const
	{
	}


	void OperatorGeneric::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
	{
	}


	string OperatorGeneric::getName() const
	{
		return "";
	}


	string OperatorGeneric::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	string OperatorGeneric::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


} // namespace acl
