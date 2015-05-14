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


#include "aclElementAssignmentSafe.h"
#include "../../aslUtilities.h"
#include "../Kernels/aclKernelConfigurationTemplates.h"

namespace acl
{


ElementAssignmentSafe::ElementAssignmentSafe(Element a1, Element a2):
	OperatorBinary(a1, a2)
{
}

string ElementAssignmentSafe::str(const KernelConfiguration & kernelConfig) const
{
	string e1s(e1->str(KERNEL_BASIC));
	string e2s(e2->str(kernelConfig));	

	if (kernelConfig.unaligned && kernelConfig.vectorWidth > 1)
	{
		return "vstore" + asl::numToStr(kernelConfig.vectorWidth) + "(" + e2s + ", 0, &" + e1s + ")";
	}
	else
	{
		return e1s + " = " + e2s;
	}
}


} // namespace acl
