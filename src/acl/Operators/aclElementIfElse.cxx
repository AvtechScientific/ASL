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


#include "aclElementIfElse.h"
#include "../aclUtilities.h"

using namespace asl;

namespace acl
{

ElementIfElse::ElementIfElse(Element condition_):
	// typeID fictitious
	ElementBase(false, 0, TYPE_INT),
	condition(condition_)
{
}


string ElementIfElse::str(const KernelConfiguration & kernelConfig) const
{
	string s("");

	// only if "if"'s body is not empty
	// create if (...) statement
	if (expressionIf.size() > 0)
	{
		s = "if (" + condition->str(kernelConfig) + ")\n\t{";

		for (unsigned int i = 0; i < expressionIf.size(); i++)
		{
			s += "\n\t\t" + expressionIf[i]->str(kernelConfig) + ";";
		}

		s += "\n\t}";

		// only if "if"'s and "else"'s body is not empty
		// create else {...} statement
		if (expressionElse.size() > 0)
		{
			s += "\n\telse\n\t{";

			for (unsigned int i = 0; i < expressionElse.size(); i++)
			{
				s += "\n\t\t" + expressionElse[i]->str(kernelConfig) + ";";
			}

			s += "\n\t}";
		}
	}
	
	return s;
}


void ElementIfElse::addBodyExpressionIf(Element expression_)
{
	if (compatible(size, queue, expression_))
	{
		size = max(size, expression_->getSize());
		if (expression_->getQueue().get() != 0)
				queue = expression_->getQueue();
		expressionIf.push_back(expression_);
	}
	else
	{
		errorMessage("ElementIfElse::addBodyExpressionIf - last added expression \
					 is incompatible with the previous ones. \
					 Either they reside on different devices or their sizes do not match: "
		             + numToStr(expression_->getSize()) + " and " + numToStr(size));
	}
}


void ElementIfElse::addBodyExpressionElse(Element expression_)
{
	if (compatible(size, queue, expression_))
	{
		size = max(size, expression_->getSize());
		if (expression_->getQueue().get() != 0)
				queue = expression_->getQueue();
		expressionElse.push_back(expression_);
	}
	else
	{
		errorMessage("ElementIfElse::addBodyExpressionElse - last added expression \
					 is incompatible with the previous ones. \
					 Either they reside on different devices or their sizes do not match: "
		             + numToStr(expression_->getSize()) + " and " + numToStr(size));
	}
}


void ElementIfElse::addToKernelSource(vector<Element> & arguments,
                                      vector<Element> & localDeclarations) const
{
	addElementToKernelSource(condition, arguments, localDeclarations);
		
	for (unsigned int i = 0; i < expressionIf.size(); i++)
	{
		addElementToKernelSource(expressionIf[i], arguments, localDeclarations);
	}

	for (unsigned int i = 0; i < expressionElse.size(); i++)
	{
		addElementToKernelSource(expressionElse[i], arguments, localDeclarations);
	}
}


void ElementIfElse::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
{
}


string ElementIfElse::getName() const
{
	return "";
}


string ElementIfElse::getTypeSignature(const KernelConfiguration & kernelConfig) const
{
	return "";
}


string ElementIfElse::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
{
	return "";
}
} // namespace acl
