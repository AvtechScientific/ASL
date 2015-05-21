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


#include "aclVectorOfElementsDef.h"
#include "aclVectorOfElementsOperations.h"
#include <aclUtilities.h>
#include <acl.h>
#include <DataTypes/aclConstant.h>
#include <../aslUtilities.h>
#include <algorithm>
#include <Kernels/aclKernel.h>

#include <acl/DataTypes/aclArray.h>
#include"../aclHardware.h"

#include "aclTypesList.h"

using asl::errorMessage;

using namespace std;

namespace acl
{
	VectorOfElementsData::VectorOfElementsData():
		std::vector<ElementData>()
	{
	}


	VectorOfElementsData::VectorOfElementsData(unsigned int n):
		std::vector<ElementData>(n)
	{
	}

	template <typename T> VectorOfElementsData::VectorOfElementsData(unsigned int n,
	                                                                 unsigned int s,
	                                                                 T,
	                                                                 CommandQueue queue):
		std::vector<ElementData>(n)
	{
		for (unsigned int i(0); i < n; ++i)
			(*this)[i] = ElementData(new Array<T>(s, queue));
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElementsData::VectorOfElementsData(unsigned int n, unsigned int s, t, CommandQueue queue);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElementsData::VectorOfElementsData(unsigned int n,
	                                                                 unsigned int s,
	                                                                 T):
		VectorOfElementsData(n,s,T(0),hardware.defaultQueue)
	{
	}
	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElementsData::VectorOfElementsData(unsigned int n, unsigned int s, t);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
		

	void VectorOfElementsData::resizeElements(unsigned int n)
	{
		if ((*this)[0].get() == 0)
			return ;
		else
			if ((*this)[0]->getSize() == n)
				return ;

		auto q((*this)[0]->getQueue());
		
		for (unsigned int i(0); i < size(); ++i)
		{
			TypeID t((*this)[i]->getTypeID());
			(*this)[i] = generateElementArray(t, n, q);
		}
	}
	

	VectorOfElements VectorOfElementsData::operator=(const VectorOfElements &a) const
	{
		return assignmentSafe(*this, a);
	}

	VectorOfElements VectorOfElementsData::operator=(const VectorOfElementsData &a) const
	{
		return assignmentSafe(*this, a);
	}
		
		

	VectorOfElements::VectorOfElements():
		std::vector<Element>()
	{
	}

		
	VectorOfElements::VectorOfElements(unsigned int n):
		std::vector<Element>(n, Element(new Constant<cl_int>(0)))
	{
	}


	VectorOfElements::VectorOfElements(const VectorOfElementsData & a):
		std::vector<Element>(a.size())
	{
		for (unsigned int i(0); i < a.size(); ++i)
			(*this)[i] = a[i];
	}
		
	
	VectorOfElements VectorOfElements::operator=(const VectorOfElements &a) const
	{
		if (size() != a.size())
			errorMessage("operator= - the sizes of two VectorOfElements are incompatible: "
			             + numToStr(size()) + " and " + numToStr(a.size()));

		VectorOfElements ve(a.size());
		for (unsigned int i(0); i < size(); ++i)
			ve[i] = elementOperators::operatorAssignment((*this)[i], a[i]);
		return ve;
	}		

	/// \todoIgal change name to: compatible()?
	bool VectorOfElements::checkCompatibility() const
	{
		unsigned int s((*this)[0]->getSize());
		CommandQueue queue((*this)[0]->getQueue());
		bool isCompatible(true);
		for (unsigned int i = 1; i < size(); ++i)
		{
			isCompatible &= compatible(s, queue, (*this)[i]);
			s = std::max(s, (*this)[i]->getSize());
			CommandQueue q((*this)[0]->getQueue());
			queue = (q.get() != 0 ? q : queue);
		}
		return isCompatible;
	}


	/// \todoIgal change name to: compatibleSizes()?
	bool VectorOfElements::checkSizesCompatibility(unsigned int n) const
	{
		bool isCompatible(true);
		for (unsigned int i = 0; i < size(); ++i)
		{
			isCompatible &= compatibleSizes(n, (*this)[i]);
		}
		return isCompatible;
	}

	unsigned int getElementsSize(const VectorOfElements & a)
	{
		if (a.size() > 0 && a[0].get() != 0)
			return a[0]->getSize();
		else
			return 0;
	}
		
	acl::TypeID getElementType(const VectorOfElements & a, unsigned int i)
	{
		return a.at(i)->getTypeID();
	}
} // namespace acl
