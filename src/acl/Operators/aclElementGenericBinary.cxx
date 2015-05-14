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


#include "aclElementGenericBinary.h"

using namespace acl;

ElementGenericBinary::ElementGenericBinary(Element a1, Element a2, const string & operation_):
	OperatorBinary(a1, a2),
	operation(operation_)
{
}


string ElementGenericBinary::str(const KernelConfiguration & kernelConfig) const
{
		return "(" + e1->str(kernelConfig) + operation + e2->str(kernelConfig) + ")";
}


ElementGenericBinaryFunction::ElementGenericBinaryFunction(Element a1, Element a2, const string & functionName_):
	OperatorBinary(a1, a2),
	functionName(functionName_)
{
}


string ElementGenericBinaryFunction::str(const KernelConfiguration & kernelConfig) const
{
	return functionName + "(" + e1->str(kernelConfig) + ", " + e2->str(kernelConfig) + ")";
}

