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


#include "aclElementSin.h"

namespace acl
{

	ElementSin::ElementSin(Element a):
		OperatorUnary(a)
	{
	}

	string ElementSin::str(const KernelConfiguration & kernelConfig) const
	{
		string es(e->str(kernelConfig));
		// check for zeros and throw them away
		if (es == "0" || es == "-0")
			return "0";
		else
			return "sin(" + es + ")";
	}

} // namespace acl
