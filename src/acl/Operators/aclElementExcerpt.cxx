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


#include "../aclUtilities.h"
#include "../../aslUtilities.h"
#include "aclElementExcerpt.h"

using namespace asl;

namespace acl
{

ElementExcerpt::ElementExcerpt(Element source_, Element filter_):
	ElementBase(source_->isWritable, filter_->getSize(), source_->getTypeID())
{
	if (onSameDevice(source_, filter_))
	{
		source = source_;
		filter = filter_;

		// filter->getQueue().get() can be 0,
		// and source->getQueue().get() should not be 0
		// (but even if it is 0 the statement bellow is correct;
		// the excerpt() operator will do nothing however in such a case).
		queue = source->getQueue();
	}
	else
	{
		errorMessage("ElementExcerpt::ElementExcerpt() - elements reside on different devices");
	}
}


string ElementExcerpt::str(const KernelConfiguration & kernelConfig) const
{
	string s = source->str(kernelConfig);
	size_t found = s.find(INDEX);

	int filterSize(filter->str(kernelConfig).size());
	while (found != string::npos)
	{
		s.replace(found, INDEX.size(), filter->str(kernelConfig));
		found = s.find(INDEX, found + filterSize);
	}

	return s;
}


void ElementExcerpt::addToKernelSource(vector<Element> & arguments,
		                               vector<Element> & localDeclarations) const
{
	addElementToKernelSource(source, arguments, localDeclarations);
	addElementToKernelSource(filter, arguments, localDeclarations);
}


void ElementExcerpt::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
{
}


string ElementExcerpt::getName() const
{
	return "";
}


string ElementExcerpt::getTypeSignature(const KernelConfiguration & kernelConfig) const
{
	return "";
}


string ElementExcerpt::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
{
	return "";
}


} // namespace acl
