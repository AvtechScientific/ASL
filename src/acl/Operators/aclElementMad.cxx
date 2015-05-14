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


#include "aclElementMad.h"

using namespace acl;


ElementMad::ElementMad(Element a1, Element a2, Element a3):
	OperatorTernary(a1, a2, a3)
{
}

string ElementMad::str(const KernelConfiguration & kernelConfig) const
{
	return "mad(" + e1->str(kernelConfig) + ", " + e2->str(kernelConfig) +", " + e3->str(kernelConfig) + ")";
}
