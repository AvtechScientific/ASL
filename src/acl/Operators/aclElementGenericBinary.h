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


#ifndef ACLELEMENTGENERICBINARY_H
#define ACLELEMENTGENERICBINARY_H

#include "../Operators/aclOperatorBinary.h"

namespace acl
{
	///\todo{Add optimization for Constant  on function level}
	///\todo{Add optimization for Singlevalued elements on function level}
	class ElementGenericBinary: public OperatorBinary
	{
		private:
			const string operation;
		public:
			ElementGenericBinary(Element a1, Element a2, const string & operation_);
			virtual string str(const KernelConfiguration & kernelConfig) const;
	};


	class ElementGenericBinaryFunction: public OperatorBinary
	{
		private:
			const string functionName;
		public:
			ElementGenericBinaryFunction(Element a1, Element a2, const string & functionName_);
			virtual string str(const KernelConfiguration & kernelConfig) const;
	};

	
} // namespace acl

#endif // ACLELEMENTGENERICBINARY_H
