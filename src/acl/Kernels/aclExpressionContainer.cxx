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


#include "aclExpressionContainer.h"
#include "../aclUtilities.h"
#include "../../aslUtilities.h"
#include <acl/aclElementBase.h>


using namespace std;
using namespace asl;

namespace acl
{

	ExpressionContainer::ExpressionContainer():
		size(0),
		queue(),
		regenerateKernelSource(true)
	{
	}


	void ExpressionContainer::addExpression(Element expression_)
	{
		if (compatible(size, queue, expression_))
		{
			size = max(size, expression_->getSize());

			if (expression_->getQueue().get() != 0)
			{
				queue = expression_->getQueue();
			}
			
			expression.push_back(expression_);
			addElementToKernelSource(expression_, arguments, localDeclarations);
			regenerateKernelSource = true;
		}
		else
		{
			errorMessage("ExpressionContainer::addExpression() - last added expression \
						 is incompatible with the previous ones. \
						 Either they reside on different devices or their sizes do not match: "
			             + numToStr(expression_->getSize()) + " and " + numToStr(size));
		}
	}

	
	void ExpressionContainer::filterDeclarations()
	{
		/// sorts all Elements in arguments
		sort(arguments.begin(), arguments.end());
		vector<Element>::iterator iter;
		/// removes consecutive duplicate Elements in arguments
		iter = unique(arguments.begin(), arguments.end());

		// drop all elements beyond the new end() (= iter) returned by unique
		arguments.resize(iter - arguments.begin());


		/// sorts all Elements in localDeclarations
		sort(localDeclarations.begin(), localDeclarations.end());

		/// removes consecutive duplicate Elements in localDeclarations
		iter = unique(localDeclarations.begin(), localDeclarations.end());
		// drop all elements beyond the new end() (= iter) returned by unique
		localDeclarations.resize(iter - localDeclarations.begin());
		
		/// \todoZeev remember duplicate Elements and create an extra local
		/// variable for them in order to reduce number of memory reads

		// Maybe use std::set to get rid of duplicates? Like this:
		// std::set<std::string>   alternatives_set (m_alternatives.begin(), m_alternatives.end());
		// std::vector<std::string> alternatives_vec (alternatives_set.begin(), alternatives_set.end());
	}

	unsigned int ExpressionContainer::getSize()
	{
		return size;
	}
	

} // namespace acl
