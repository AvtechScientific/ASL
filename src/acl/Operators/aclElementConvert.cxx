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


#include "aclElementConvert.h"
#include <aslUtilities.h>
#include "../aclHardware.h"

namespace acl
{

	ElementConvert::ElementConvert(Element a, TypeID type, bool strong):
		OperatorUnary(a,type),
		strongForm(strong)
	{
	}
		
	string ElementConvert::str(const KernelConfiguration & kernelConfig) const
	{
		string es(e->str(kernelConfig));
		string vectorW(kernelConfig.vectorWidth>1 ? 
		               asl::numToStr(kernelConfig.vectorWidth) :
			           "");
		string function;
		
		if (strongForm)
			function="convert_"+TYPE[typeID]+vectorW;
		else
			function="("+TYPE[typeID]+vectorW+")";
		return function+"(" + es + ")";
	}


} // namespace acl
