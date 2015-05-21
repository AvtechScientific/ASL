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


#include "aclElementGenericUnary.h"

using namespace acl;

ElementGenericUnary::ElementGenericUnary(Element a, 
                                         const string & operation_, 
                                         bool outside):
	OperatorUnary(a),
	operation(operation_),
	outside(outside)
{
}

string ElementGenericUnary::str(const KernelConfiguration & kernelConfig) const
{
	string res;
	if (outside)
		res = operation + "(" + e->str(kernelConfig) + ")";
	else 
		res =  "(" + operation + e->str(kernelConfig) + ")";
	return res;
}

ElementGenericUnarySIMD::ElementGenericUnarySIMD(Element a, 
                                         const string & operation_, 
                                         bool outside):
	OperatorUnary(a),
	operation(operation_),
	outside(outside)
{
}

string ElementGenericUnarySIMD::str(const KernelConfiguration & kernelConfig) const
{
	string res;
	if(kernelConfig.vectorWidth > 1)
	{
		if (outside)
			res = operation + "(" + e->str(kernelConfig) + ")";
		else 
			res =  "(" + operation + e->str(kernelConfig) + ")";
	}
	else 
		res = e->str(kernelConfig);	
	return res;
}
