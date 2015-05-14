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


#include "aclElementParser.h"
#include "../aclUtilities.h"
#include "../../aslUtilities.h"

using namespace asl;

namespace acl
{

ElementParser::ElementParser():
	// set typeID fictitious first, change with the first call of addElementNamePair
	ElementBase(false, 0, TYPE_INT)
{
}


void ElementParser::addElementNamePair(Element element, string name)
{
	if (compatible(size, queue, element))
	{
		size = max(size, element->getSize());
		typeID = element->getTypeID();
		if (element->getQueue().get() != 0)
				queue = element->getQueue();
		elementNamePairs.push_back(make_pair(element, name));
	}
	else
	{
		errorMessage("ElementParser::addElementNamePair - last added expression \
					 is incompatible with the previous ones. \
					 Either they reside on different devices or their sizes do not match: "
		             + numToStr(element->getSize()) + " and " + numToStr(size));
	}
}


void ElementParser::setStatement(string statement_)
{
	statement = statement_;
}


string ElementParser::str(const KernelConfiguration & kernelConfig) const
{
	string s = statement;

	for (unsigned int i = 0; i < elementNamePairs.size(); ++i)
	{
		size_t found = s.find(elementNamePairs[i].second);

		int nameSize(elementNamePairs[i].second.size());
		string fieldStr(elementNamePairs[i].first->str(kernelConfig));
		while (found != string::npos)
		{
			s.replace(found, nameSize, fieldStr);
			found = s.find(elementNamePairs[i].second, found + fieldStr.size());
		}
	}

	return s;
}

/// \todoZeev Add to kernel only Elements that were found in the statement!!!
void ElementParser::addToKernelSource(vector<Element> & arguments,
		                              vector<Element> & localDeclarations) const
{
	if (statement == "")
		errorMessage("ElementParser::addToKernelSource() - statement is an empty string");
	
	for (unsigned int i = 0; i < elementNamePairs.size(); ++i)
	{
		addElementToKernelSource(elementNamePairs[i].first, arguments, localDeclarations);
	}
}


void ElementParser::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
{
}


string ElementParser::getName() const
{
	return "";
}


string ElementParser::getTypeSignature(const KernelConfiguration & kernelConfig) const
{
	return "";
}


string ElementParser::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
{
	return "";
}
} // namespace acl
