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


#ifndef ACLKERNELCONFIGURATION_H
#define ACLKERNELCONFIGURATION_H

#include <string>
#include <vector>

namespace acl
{

	/// ACL Kernel configuration class 	\ingroup KernelGen
	class KernelConfiguration
	{
		public:
			explicit KernelConfiguration(bool simd = false,
			                             bool unaligned_ = false,
			                             bool local_ = false);
			KernelConfiguration(const KernelConfiguration & kernelConfig_);
			bool operator==(const KernelConfiguration & a) const;

			unsigned int vectorWidth;
			bool unaligned;
			bool local;
			std::vector<std::string> extensions;
	};


} // acl namespace

#endif // ACLKERNELCONFIGURATION_H
