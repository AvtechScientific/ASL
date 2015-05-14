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


#include "aclElementFor.h"
#include "../aclUtilities.h"
#include "../../aslUtilities.h"

using namespace asl;

namespace acl
{

ElementFor::ElementFor(Element initialization_, Element condition_, Element increase_):
	// typeID fictitious
	ElementBase(false, 0, TYPE_INT),
	initialization(initialization_),
	condition(condition_),
	increase(increase_)
{
	if (compatible(size, queue, initialization_))
	{
		size = max(size, initialization_->getSize());
		if (initialization_->getQueue().get() != 0)
			queue = initialization_->getQueue();
	}
	if (compatible(size, queue, condition_))
	{
		size = max(size, condition_->getSize());
		if (condition_->getQueue().get() != 0)
			queue = condition_->getQueue();
	}
	if (compatible(size, queue, increase_))
	{
		size = max(size, increase_->getSize());
		if (increase_->getQueue().get() != 0)
			queue = increase_->getQueue();
	}
}


string ElementFor::str(const KernelConfiguration & kernelConfig) const
{
	string s("");

	// only if "for"'s body is not empty
	// create for (...) loop
	if (expression.size() > 0)
	{
		s = "for (" + initialization->str(kernelConfig) + "; "
			+ condition->str(kernelConfig) + "; "
			+ increase->str(kernelConfig) + ")\n\t{";

		for (unsigned int i = 0; i < expression.size(); i++)
		{
			s += "\n\t\t" + expression[i]->str(kernelConfig) + ";";
		}

		s += "\n\t}";
	}
	
	return s;
}


void ElementFor::addBodyExpression(Element expression_)
{
	if (compatible(size, queue, expression_))
	{
		size = max(size, expression_->getSize());
		if (expression_->getQueue().get() != 0)
				queue = expression_->getQueue();
		expression.push_back(expression_);
	}
	else
	{
		errorMessage("ElementFor::addBodyExpression - last added expression \
					 is incompatible with the previous ones. \
					 Either they reside on different devices or their sizes do not match: "
		             + numToStr(expression_->getSize()) + " and " + numToStr(size));
	}
}


void ElementFor::addToKernelSource(vector<Element> & arguments,
		                           vector<Element> & localDeclarations) const
{
	addElementToKernelSource(initialization, arguments, localDeclarations);
	addElementToKernelSource(condition, arguments, localDeclarations);
	addElementToKernelSource(increase, arguments, localDeclarations);
	
	for (unsigned int i = 0; i < expression.size(); i++)
	{
		addElementToKernelSource(expression[i], arguments, localDeclarations);
	}
}


void ElementFor::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
{
}


string ElementFor::getName() const
{
	return "";
}


string ElementFor::getTypeSignature(const KernelConfiguration & kernelConfig) const
{
	return "";
}


string ElementFor::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
{
	return "";
}
} // namespace acl
