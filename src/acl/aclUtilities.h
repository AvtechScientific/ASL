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


#ifndef ACLUTILITIES_H
#define ACLUTILITIES_H

#include "aslUtilities.h"
#include "Kernels/aclKernelConfiguration.h"
//#include "aclHardware.h"
#include "aclTypes.h"

namespace cl
{
	class CommandQueue;
}

namespace acl
{
	extern const std::string INDEX;
	extern const KernelConfiguration KERNEL_BASIC;

	class ElementBase;
	typedef std::shared_ptr<ElementBase> Element;
	typedef std::shared_ptr<cl::CommandQueue> CommandQueue;
	
	/// determines whether two elements are compatible
	/// i.e. have compatible sizes and reside on the same device
	bool compatible(unsigned int size1, CommandQueue queue1,
	                unsigned int size2, CommandQueue queue2);
	bool compatible(unsigned int size, CommandQueue queue, Element e);
	bool compatible(Element e1, Element e2);

	/// if \p s and size of \p e are the same or at least one of them is 0
	inline bool compatibleSizes(unsigned int s1, unsigned int s2);
	/// if both elements have the same size or size of at least one of them is 0
	bool compatibleSizes(Element e1, Element e2);
	/// if \p s and size of \p e are the same or at least one of them is 0
	bool compatibleSizes(unsigned int s, Element e);

	/// Adds padding in bytes based on the device's alignment
	unsigned int paddingBytes(unsigned int size,
	                          unsigned int typeSize,
	                          CommandQueue queue);

	/// Adds padding in elements based on vector width
	unsigned int paddingElements(unsigned int size,
	                             const KernelConfiguration & kernelConfig);
	
	/// checks whether both elements reside on the same device
	bool onSameDevice(CommandQueue queue1, CommandQueue queue2);
	bool onSameDevice(CommandQueue queue, Element e);
	bool onSameDevice(Element e1, Element e2);

	bool isDeclarable(Element e);

	bool isArgument(Element e);

	/// adds \p e either to \p arguments or to \p localDeclarations
	void addElementToKernelSource(Element e,
	                              std::vector<Element> & arguments,
			                      std::vector<Element> & localDeclarations);
	
	template <typename T> const std::string& typeToStr();
	template <typename T> inline const std::string typeToStr(unsigned int i);

	template <typename T> inline constexpr const TypeID typeToTypeID();
	
// ------------------------Implementation ---------------

	inline bool compatibleSizes(unsigned int s1, unsigned int s2)
	{
		return (s1 == s2 || (s1 * s2) == 0);
	}

	template <> inline constexpr const TypeID typeToTypeID<cl_double>()
	{
		return TYPE_DOUBLE;
	}

	template <> inline constexpr const TypeID typeToTypeID<cl_float>()
	{
		return TYPE_FLOAT;
	}

	template <> inline constexpr const TypeID typeToTypeID<cl_int>()
	{
		return TYPE_INT;
	}

	template <> inline constexpr const TypeID typeToTypeID<cl_uint>()
	{
		return TYPE_UINT;
	}

	template <> inline constexpr const TypeID typeToTypeID<cl_long>()
	{
		return TYPE_LONG;
	}
		
	template <typename T> inline const std::string typeToStr(unsigned int i)
	{
		if (i == 1)
			return typeToStr<T>();
		else
			return typeToStr<T>() + asl::numToStr(i);		
	}
	
	
} // namespace acl
#endif // ACLUTILITIES_H