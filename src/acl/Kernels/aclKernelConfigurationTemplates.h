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


#ifndef ACLKERNELCONFIGURATIONTEMPLATES_H
#define ACLKERNELCONFIGURATIONTEMPLATES_H

#include "aclKernelConfiguration.h"

namespace acl
{

	// GLOBALS


	/// \ingroup KernelGen
	extern const KernelConfiguration KERNEL_BASIC;
	/// \ingroup KernelGen
	extern const KernelConfiguration KERNEL_SIMD;
	/// \ingroup KernelGen
	extern const KernelConfiguration KERNEL_SIMDUA;
	

} // acl namespace

#endif // ACLKERNELCONFIGURATIONTEMPLATES_H
