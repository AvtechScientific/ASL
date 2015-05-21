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


#include "aclElementDivision.h"
#include "../../aslUtilities.h"

namespace acl
{


ElementDivision::ElementDivision(Element a1, Element a2):
	OperatorBinary(a1, a2)
{
}


string ElementDivision::str(const KernelConfiguration & kernelConfig) const
{
	string e1s(e1->str(kernelConfig));
	string e2s(e2->str(kernelConfig));

	// check for zeros and throw them away
	if ((e1s != "0") && (e2s != "0") && (e1s != "-0") && (e2s != "-0"))
		return "(" + e1s + "/" + e2s + ")";
	else
		if (((e1s == "0") || (e1s == "-0")) && (e2s != "0") && (e2s != "-0"))
		    return "0";
		else
			asl::errorMessage("Division by zero");
    return "";
}

} // namespace acl
