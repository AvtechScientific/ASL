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


#ifndef ACLELEMENTCONVERT_H
#define ACLELEMENTCONVERT_H

#include "../Operators/aclOperatorUnary.h"


namespace acl
{
	///realizes \p convert_ functionality
	/**
		\ingroup SpecialPurposeFunctions
		The element generates code corresponding to convert function with 
		 an apropriate vector index.

		The \p strong parameter defines the form of the operator: 
		 convert_ in case of \p true or (type) in case of false 
	*/
	class ElementConvert: public OperatorUnary
	{
		private: 
			bool strongForm;
		public:
			ElementConvert(Element a, TypeID type, bool strong=true);
			virtual string str(const KernelConfiguration & kernelConfig) const;
	};

} // namespace acl

#endif // ACLELEMENTCONVERT_H
