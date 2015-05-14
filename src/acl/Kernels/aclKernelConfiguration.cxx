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


#include "aclKernelConfiguration.h"


using namespace std;

namespace acl
{

KernelConfiguration::KernelConfiguration(bool simd,
                                         bool unaligned_,
                                         bool local_):
	// sets vectorWidth to 2 for now, later it will be automatically
	// detected and overwritten by the Kernel
	vectorWidth(simd ? 2 : 1),
	unaligned(unaligned_),
	local(local_)

{
}


KernelConfiguration::KernelConfiguration(const KernelConfiguration & kernelConfig_):
	vectorWidth(kernelConfig_.vectorWidth),
	unaligned(kernelConfig_.unaligned),
	local(kernelConfig_.local)

{
	extensions = kernelConfig_.extensions;
}


bool KernelConfiguration::operator==(const KernelConfiguration & a) const
{
	bool equal;

	equal = (vectorWidth == a.vectorWidth)
			&& (unaligned == a.unaligned)
			&& (local == a.local)
			&& (extensions == a.extensions);

	return equal;

}



} // asl namespace
