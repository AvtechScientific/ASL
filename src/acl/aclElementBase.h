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


#ifndef ACLELEMENTBASE_H
#define ACLELEMENTBASE_H

#include "aclStdIncludes.h"
//#include "aclHardware.h"
#include "Kernels/aclKernelConfiguration.h"
#include <memory>
#include "aclTypes.h"

namespace cl
{
	class CommandQueue;
	class Kernel;
}

namespace acl
{
	typedef std::shared_ptr<cl::CommandQueue> CommandQueue;

	// GLOBALS	
	extern const std::string INDEX;
	extern const KernelConfiguration KERNEL_BASIC;
	
	using namespace std;
	
	class ElementBase
	{
		protected:
			unsigned int size;
			TypeID typeID;
			CommandQueue queue;
			ElementBase(bool isWritable_, unsigned int size_, TypeID typeID_);
		public:
			const bool isWritable;
			virtual string str(const KernelConfiguration & kernelConfig = KERNEL_BASIC) const = 0;
			virtual string getName() const = 0;
			virtual string getAddressSpaceQualifier() const;
			unsigned int getSize() const;
			CommandQueue getQueue() const;
			TypeID getTypeID() const;
			virtual string getTypeSignature(const KernelConfiguration & kernelConfig = KERNEL_BASIC) const = 0;
			virtual string getLocalDeclaration(const KernelConfiguration & kernelConfig = KERNEL_BASIC) const = 0;
			/// Adds ElementBase to the kernel source either as an argument or as a local declaration
			virtual void addToKernelSource(vector<shared_ptr<ElementBase> > & arguments,
			                               vector<shared_ptr<ElementBase> > & localDeclarations) const = 0;

			virtual void setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const = 0;
			virtual ~ElementBase();
	};

	typedef std::shared_ptr<ElementBase> Element;
	
} // namespace acl
#endif // ACLELEMENTBASE_H
