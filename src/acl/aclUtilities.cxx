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


#include "aclUtilities.h"
#include "Kernels/aclKernelConfigurationTemplates.h"
#include "acl/aclElementBase.h"

#include "aclTypesList.h"
#include "aclHardware.h"

using namespace std;

namespace acl
{

	bool compatible(unsigned int size1,
	                CommandQueue queue1,
	                unsigned int size2,
	                CommandQueue queue2)
	{
		return ((compatibleSizes(size1, size2)) && (onSameDevice(queue1, queue2)));
	}


	bool compatible(unsigned int size, CommandQueue queue, Element e)
	{
		return compatible(size, queue, e->getSize(), e->getQueue());
	}


	bool compatible(Element e1, Element e2)
	{
		return compatible(e1->getSize(), e1->getQueue(), e2->getSize(), e2->getQueue());
	}

	
	bool compatibleSizes(Element e1, Element e2)
	{
		return compatibleSizes(e1->getSize(), e2->getSize());
	}


	bool compatibleSizes(unsigned int s, Element e)
	{
		return compatibleSizes(s, e->getSize());
	}


	bool onSameDevice(CommandQueue queue1, CommandQueue queue2)
	{
		return ((queue1 == queue2) || (queue1.get() == 0) || (queue2.get() == 0));
	}

	bool onSameDevice(CommandQueue queue, Element e)
	{
		return onSameDevice(queue, e->getQueue());
	}

	
	bool onSameDevice(Element e1, Element e2)
	{
		return onSameDevice(e1->getQueue(), e2->getQueue());
	}


	bool isDeclarable(Element e)
	{
		return ((e->getTypeSignature(KERNEL_BASIC) != "") || (e->getLocalDeclaration(KERNEL_BASIC) != ""));
	}


	bool isArgument(Element e)
	{
		return (e->getTypeSignature(KERNEL_BASIC) != "");
	}

	
	void addElementToKernelSource(Element e,
	                              vector<Element> & arguments,
			                      vector<Element> & localDeclarations)
	{
		// adds expression_ either to arguments or to localDeclarations
		// in a case that expression_ isDeclarable
		if (isDeclarable(e))
		{
			if (isArgument(e))
			{
				arguments.push_back(e);
			}
			else
			{
				localDeclarations.push_back(e);
			}
		}
		else
		{
			// adds expressions' elements either to arguments or to localDeclarations;
			// multiply occurences of the same element possible.
			// needs to be filtered
			e->addToKernelSource(arguments, localDeclarations);
		}
	}


	unsigned int paddingBytes(unsigned int size, unsigned int typeSize, CommandQueue queue)
	{
		// second modulo is added in order to make padding = 0 in the case
		// that (size * typeSize) is divisible by getAlignment(queue)
		return (getAlignment(queue) - ((size * typeSize)  % getAlignment(queue))) % getAlignment(queue);
	}


	unsigned int paddingElements(unsigned int size, const KernelConfiguration & kernelConfig)
	{
		// second modulo is added in order to make padding = 0 in the case
		// that size is divisible by kernelConfig.vectorWidth
		return (kernelConfig.vectorWidth - (size  % kernelConfig.vectorWidth)) % kernelConfig.vectorWidth;
	}

	template <typename T> const std::string& typeToStr()
	{
		return TYPE[typeToTypeID<T>()];
	}

	#define BOOST_TT_rep_expression(r, data, T) \
	template const std::string& typeToStr<T>();
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	
} // namespace acl