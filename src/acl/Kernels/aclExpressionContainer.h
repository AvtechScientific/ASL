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


#ifndef ACLEXPRESSIONCONTAINER_H
#define ACLEXPRESSIONCONTAINER_H

#include <vector>
#include <memory>
//#include <CL/cl.hpp>
// Supply "cl.hpp" with ASL, since it is not present in OpenCL 2.0
// Remove the file after switching to OpenCL 2.1
#include "acl/cl.hpp"


namespace acl
{

	typedef std::shared_ptr<cl::CommandQueue> CommandQueue;

	class ElementBase;
	typedef std::shared_ptr<ElementBase> Element;

	class ExpressionContainer
	{
		protected:
			unsigned int size; //< initialized with 0, is overidden by first addExpression()
			CommandQueue queue;
			bool regenerateKernelSource;
		public:	
			std::vector<Element> expression;
		protected:
			std::vector<Element> arguments;
			std::vector<Element> localDeclarations;
			void filterDeclarations();
		public:
			ExpressionContainer();
			void addExpression(Element expression_);		
			unsigned int getSize();
			inline const CommandQueue getQueue () const; 
	};


//------------------------------ Imlementation --------------------

	inline const CommandQueue ExpressionContainer::getQueue () const 
	{
		return queue;
	}
	
} // namespace acl

#endif // ACLEXPRESSIONCONTAINER_H
