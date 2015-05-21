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


#include "aclElementSyncCopy.h"
#include "../../aslUtilities.h"
#include "../aclUtilities.h"
#include "../aclHardware.h"

using namespace asl;

namespace acl
{

///\todo{Think of cutting one copy operation into 3: frontRemainder, middlePart (SIMD chunks) and backRemainder}
	
unsigned int ElementSyncCopy::syncCopyNum = 0;


ElementSyncCopy::ElementSyncCopy(Element source_,
                                 Element destination_,
                                 Element srcOffset,
                                 Element dstOffset,
                                 Element length_):
	// typeID fictitious
	ElementBase(false, 0, TYPE_INT),
	source(source_),
	destination(destination_),
	length(length_)
{
	if(source->getQueue().get() == 0 && destination->getQueue().get() == 0)
		errorMessage("ElementSyncCopy: source and destionation have no queue");
	queue = source->getQueue().get() != 0 ? source->getQueue() : destination->getQueue();
	
	kernelSourceStr = "event_t event_" + numToStr(syncCopyNum) + " = (event_t)0;\n";
	kernelSourceStr += "\tevent_" + numToStr(syncCopyNum) + " = async_work_group_copy(&((" + destination->getAddressSpaceQualifier() + " " + TYPE[destination->getTypeID()] + " *)" + destination->getName() + ")[" + 
						dstOffset->str() + "], &((" + source->getAddressSpaceQualifier() + " " + TYPE[source->getTypeID()] + " *)" + source->getName() + ")[" +
						srcOffset->str() + "], " +	length->str() +
						", (event_t)0);\n";
	kernelSourceStr += "\twait_group_events (1, &event_" + numToStr(syncCopyNum) + ")";

	++syncCopyNum;
}


string ElementSyncCopy::str(const KernelConfiguration & kernelConfig) const
{
	return kernelSourceStr;
}


void ElementSyncCopy::addToKernelSource(vector<Element> & arguments,
                                        vector<Element> & localDeclarations) const
{
	if (isDeclarable(source))
	{
		if (isArgument(source))
		{
			arguments.push_back(source);
		}
		else
		{
			localDeclarations.push_back(source);
		}
	}
	else
	{
		source->addToKernelSource(arguments, localDeclarations);
	}
		
	if (isDeclarable(destination))
	{
		if (isArgument(destination))
		{
			arguments.push_back(destination);
		}
		else
		{
			localDeclarations.push_back(destination);
		}
	}
	else
	{
		destination->addToKernelSource(arguments, localDeclarations);
	}
}


void ElementSyncCopy::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
{
}


string ElementSyncCopy::getName() const
{
	return "";
}


string ElementSyncCopy::getTypeSignature(const KernelConfiguration & kernelConfig) const
{
	return "";
}


string ElementSyncCopy::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
{
	return "";
}

} // namespace acl
