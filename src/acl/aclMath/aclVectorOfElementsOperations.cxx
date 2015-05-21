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


#include "aclVectorOfElementsOperations.h"
#include "aclUtilities.h"
#include "acl.h"
#include "DataTypes/aclConstant.h"
#include "../aslUtilities.h"
#include <algorithm>
#include "Kernels/aclKernel.h"
#include "Kernels/aclKernelConfigurationTemplates.h"

#include <acl/DataTypes/aclArray.h>

#include <acl/aclGenerators.h>

using asl::errorMessage;

using namespace std;

namespace acl
{

	VectorOfElementsData clone(VectorOfElementsData a)
	{
		if (a.size() == 0)
			return VectorOfElementsData(0u);

		VectorOfElementsData v(a.size());
		auto queue(a[0]->getQueue());
		
		for (unsigned int i(0); i < a.size(); ++i)
			v[i] = generateElementArray(a[i]->getTypeID(), a[i]->getSize(), queue);
		return v;
	}

	VectorOfElementsData clone(VectorOfElementsData a, unsigned int n)
	{
		if (a.size() < n)
			errorMessage("clone: number n is larger than size of the input vector");
		
		if (a.size() == 0)
			return VectorOfElementsData(0u);

		VectorOfElementsData v(n);
		auto queue(a[0]->getQueue());
		
		for (unsigned int i(0); i < n; ++i)
			v[i] = generateElementArray(a[i]->getTypeID(), a[i]->getSize(), queue);
		return v;
	}
	


	VectorOfElements assignmentSafe(const VectorOfElements & a,
	                                const VectorOfElements &b)
	{
		if (a.size() != b.size())
			errorMessage("assignmentSafe - the sizes of two VectorOfElements are incompatible: "
			             + asl::numToStr(a.size()) + " and " + asl::numToStr(b.size()));

		VectorOfElements ve(a.size());
		for (unsigned int i(0); i < a.size(); ++i)
			ve[i] = elementOperators::operatorAssignmentSafe(a[i], b[i]);

		return ve;
	}

	void copy(const vector<Element> & source,
	          VectorOfElements & destination)
	{
		destination.resize(source.size());
		for (unsigned int i(0); i < destination.size(); ++i)
			destination[i] = source[i];
	}

	void copy(const vector<ElementData> & source,
	          VectorOfElements & destination)
	{
		destination.resize(source.size());
		for (unsigned int i(0); i < destination.size(); ++i)
			destination[i] = source[i];
	}

	
	void copy(const vector<Element> & source,
	          VectorOfElements & destination, unsigned int start, unsigned int end)
	{
		if (source.size() <= end )
			errorMessage("copy: attempt to copy besides the vector range");
		
		destination.resize(1+end-start);
		for (unsigned int i(start); i <= end; ++i)
			destination[i] = source[i];
	}

	void copy(const VectorOfElementsData & source,
	          VectorOfElementsData & destination)
	{
		destination.resize(source.size());
		for (unsigned int i(0); i < destination.size(); ++i)
			destination[i] = source[i];
	}


	VectorOfElements subVE(const VectorOfElements & source,
	                       unsigned int start, unsigned int end)
	{
		if (source.size() <= end )
			errorMessage("subVE: attempt to copy besides the vector range");
		
		VectorOfElements destination(1 + end - start);
		for (unsigned int i(start); i <= end; ++i)
			destination[i - start] = source[i];
		return destination;
	}

	VectorOfElements subVE(const VectorOfElements & source,
	                       unsigned int i)
	{
		if (source.size() <= i )
			errorMessage("subVE: attempt to copy besides the vector range");
		
		VectorOfElements destination(1);
			destination[0] = source[i];
		return destination;
	}

	VectorOfElements subVE(const VectorOfElements & source,
	                       const vector<unsigned int> & iList)
	{
		VectorOfElements destination(iList.size());
		for (unsigned int i(0); i < iList.size(); ++i)
		{
			if (source.size() <= iList[i] )
				errorMessage("subVE: attempt to copy besides the vector range");
			destination[i] = source[iList[i]];
		}
		return destination;
	}
	
		
	VectorOfElementsData subVE(const VectorOfElementsData & source,
	                           unsigned int start, unsigned int end)
	{
		if (source.size() <= end )
			errorMessage("subVE: attempt to copy besides the vector range");
		
		VectorOfElementsData destination(1 + end - start);
		for (unsigned int i(start); i <= end; ++i)
			destination[i - start] = source[i];
		return destination;
	}

	VectorOfElements cat(const VectorOfElements * a, 
	                     unsigned int n)
	{
		unsigned int s(0);
		for(unsigned int j(0); j < n; ++j)
			s+=a[j].size();

		VectorOfElements c(s);
		s=0;
		for(unsigned int i(0); i < n; ++i)
			for(unsigned int j(0); j < a[i].size(); ++j)
			{
				c[s] = a[i][j];
				s++;
			}
		return c;
	}
	

	VectorOfElements catN(const VectorOfElements & a, 
	                      unsigned int n)
	{
		VectorOfElements c(a.size()*n);
		for(unsigned int j(0); j < a.size(); ++j)
			for(unsigned int i(0); i < n; ++i)
				c[j+i*a.size()] = a[j];		
		return c;
	}

	void initData(VectorOfElements a,
	              VectorOfElements initializationValue)
	{
		initData(a, initializationValue, KERNEL_SIMD);
	}

	
	void initData(VectorOfElements a,
	              VectorOfElements initializationValue,
	              const KernelConfiguration & kernelConfig)
	{
		Kernel k(kernelConfig);
		k << acl::assignmentSafe(a, initializationValue);
		k.setup();
		k.compute();
	}


	/// \todoIgal do they have to be on the same device?
	void swapBuffers(const VectorOfElementsData & a, const VectorOfElementsData & b)
	{ 
		if (a.size() == 0)
			errorMessage("swapBuffers - first VectorsOfEllementsData has zero size");		
		if (b.size() == 0)
			errorMessage("swapBuffers - second VectorsOfEllementsData has zero size");		
		if (a.size() != b.size())
			errorMessage("swapBuffers - twoVectorOfElementsData have different sizes");

		for (unsigned int i(0); i < a.size(); ++i)
			swapBuffers(*a[i], *b[i]);
	}
		

	VectorOfElements  operator-(const VectorOfElements & a)
	{ 
		VectorOfElements c(a.size());
		using namespace elementOperators;			
		for (unsigned int i(0); i < c.size(); ++i)
		{
			c[i]= - a[i];
		}
		return c;
	}


	VectorOfElements  operator+=(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		if (a.size() != b.size())
			errorMessage("operator+= - two VectorOfElements have different sizes");		

		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i]=( a[i] += b[i]);
		return c;
	}


	VectorOfElements  operator-=(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		VectorOfElements c(a.size());
		/// \todoIgal on same device? and all the rest of operators...
		if (a.size() == b.size())
		{
			using namespace elementOperators;			
			for (unsigned int i(0); i < c.size(); ++i)
				c[i] = (a[i] -= b[i]);
		}
		else
			errorMessage("operator-= - two VectorOfElements have different sizes");		
		return c;
	}


	VectorOfElements  operator*=(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		VectorOfElements c(a.size());
		if (b.size() == 1)
		{	
			using namespace elementOperators;
			for (unsigned int i(0); i < a.size(); ++i)
				c[i] = (a[i] *= b[0]);
		}
		else
			errorMessage("operator*= - the second VectorOfElements has more than 1 element");		
		return c;
	}


	VectorOfElements  operator/=(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		VectorOfElements c(a.size());
		if (b.size() == 1)
		{
			using namespace elementOperators;
			for (unsigned int i(0); i < a.size(); ++i)
				c[i] = (a[i] /= b[0]);
		}
		else
			errorMessage("operator/= - the second VectorOfElements has more than 1 element");		
		return c;
	}


	VectorOfElements  operator+(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		if (a.size() != b.size())
			errorMessage("operator+ - two VectorOfElements have different sizes");		

		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = (a[i] + b[i]);
		return c;
	}


	VectorOfElements  operator-(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		if (a.size() != b.size())
			errorMessage("operator- - two VectorOfElements have different sizes:"+ 
			             asl::numToStr(a.size()) + ", " +
			             asl::numToStr(b.size()));
		
		VectorOfElements c(a.size());
		using namespace elementOperators;			
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = ( a[i] - b[i]);					
		return c;
	}


	VectorOfElements  operator*(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		VectorOfElements c(1);
		if (a.size() == b.size())
		{
			using namespace elementOperators;
			c[0] = a[0] * b[0];
			for (unsigned int i(1); i < a.size(); ++i)
				c[0] = c[0] + a[i] * b[i];
			return c;
		}
		if (a.size() == 1)
		{
			c.resize(b.size());	
			using namespace elementOperators;
			for (unsigned int i(0); i < c.size(); ++i)
				c[i] = a[0] * b[i];
			return c;
		}
		if (b.size() == 1)
		{
			c.resize(a.size());	
			using namespace elementOperators;
			for (unsigned int i(0); i < c.size(); ++i)
				c[i] = a[i] * b[0];
			return c;
		}
		
		errorMessage("operator* - two VectorOfElements have different sizes or nor of one has size 1");		
		return c;
	}


	VectorOfElements  operator/(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		if (b.size() != 1u)
			errorMessage("operator/ - size of second VectorOfElements is not 1");		

		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i]= a[i] / b[0];
		return c;
	}

	VectorOfElements  operator%(const VectorOfElements & a, 
	                            const VectorOfElements & b)
	{ 
		if (b.size() != 1u)
			errorMessage("operator% - size of second VectorOfElements is not 1");		

		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i]= a[i] % b[0];
		return c;
	}

	VectorOfElements operator==(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		VectorOfElements c(1);
		if (a.size() == b.size())
		{
			using namespace elementOperators;
			c[0] = isEqual(a[0], b[0]);
			for (unsigned int i(1); i < c.size(); ++i)
				c[0] = (c[0] && isEqual(a[i], b[i]));
			return c;
		}
		else
			errorMessage("operator== - two VectorOfElements have different sizes");		
		return c;
	}


	VectorOfElements  operator!=(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		VectorOfElements c(1);
		if (a.size() == b.size())
		{
			using namespace elementOperators;
			c[0] = isNotEqual(a[0], b[0]);
			for (unsigned int i(1); i < c.size(); ++i)
				c[0] = (c[0] || isNotEqual(a[i], b[i]));
			return c;
		}		
		errorMessage("operator!=  - two VectorOfElements have different sizes");		
		return c;
	}


	VectorOfElements  operator>(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		VectorOfElements c(1);
		if (a.size() == b.size())
		{
			using namespace elementOperators;
			c[0] = (a[0] > b[0]);
			for (unsigned int i(1); i < a.size(); ++i)
				c[0] = (c[0] && a[i] > b[i]);
			return c;
		}		
		errorMessage("operator> - two VectorOfElements have different sizes");		
		return c;
	}

	VectorOfElements  operator&&(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		if (a.size() != b.size())
			errorMessage("operator&& - two VectorOfElements have different sizes");		
			
		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = a[i] && b[i];
		return c;
	}

	VectorOfElements  operator||(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		if (a.size() != b.size())
			errorMessage("operator|| - two VectorOfElements have different sizes");		
			
		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = a[i] || b[i];
		return c;
	}
	
	VectorOfElements  operator<(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		if (a.size() != b.size())
			errorMessage("operator< - two VectorOfElements have different sizes");			

		VectorOfElements c(1);		
		using namespace elementOperators;
		c[0] = (a[0] < b[0]);
		for (unsigned int i(1); i < a.size(); ++i)
			c[0] = (c[0] && a[i] < b[i]);
		return c;
	}

	VectorOfElements  operator<=(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		if (a.size() != b.size())
			errorMessage("operator<= - two VectorOfElements have different sizes");			

		VectorOfElements c(1);		
		using namespace elementOperators;
		c[0] = (a[0] <= b[0]);
		for (unsigned int i(1); i < a.size(); ++i)
			c[0] = (c[0] && a[i] <= b[i]);
		return c;
	}

	VectorOfElements  operator>=(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		if (a.size() != b.size())
			errorMessage("operator>= - two VectorOfElements have different sizes");			

		VectorOfElements c(1);		
		using namespace elementOperators;
		c[0] = (a[0] >= b[0]);
		for (unsigned int i(1); i < a.size(); ++i)
			c[0] = (c[0] && a[i] >= b[i]);
		return c;
	}
	
	VectorOfElements  operator!(const VectorOfElements & a)
	{ 
		VectorOfElements c(a.size());
		using namespace elementOperators;			
		for (unsigned int i(0); i < c.size(); ++i)
			c[i]= !a[i];
		return c;
	}
	
	VectorOfElements crossProduct(const VectorOfElements & a, const VectorOfElements & b)
	{
		if (a.size() != b.size())
			errorMessage("crossProduct - two VectorOfElements have different sizes");		
		if (a.size() == 1)
			errorMessage("crossProduct - is undefined for size 1 of VectorOfElements");		
		if (a.size() > 3)
			errorMessage("crossProduct - is undefined for size more than 3 of VectorOfElements");		

		VectorOfElements c(a.size()==2 ? 1 :3);
		using namespace elementOperators;
		
		if (a.size() == 2)
			c[0] = a[0] * b[1] - a[1] * b[0];
		else
		{
			c[0] = a[1] * b[2] - a[2] * b[1];
			c[1] = a[2] * b[0] - a[0] * b[2];
			c[2] = a[0] * b[1] - a[1] * b[0];	
		}
		return c;
	}
	

	VectorOfElements productOfElements(const VectorOfElements & a, const VectorOfElements & b)
	{
		if (a.size() != b.size())
			errorMessage("productOfElements - two VectorOfElements have different sizes");		

		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = a[i] * b[i];
		return c;
	}

	VectorOfElements productOfElements(const VectorOfElements & a)
	{
		VectorOfElements c(1);
		using namespace elementOperators;
		c[0]=a[0];
		for (unsigned int i(1); i < a.size(); ++i)
			c[0] = c[0] * a[i];
		return c;
	}

	VectorOfElements divisionOfElements(const VectorOfElements & a, const VectorOfElements & b)
	{
		if (a.size() != b.size())
			errorMessage("divisionOfElements - two VectorOfElements have different sizes");		

		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = a[i] / b[i];
		return c;
	}

	VectorOfElements min(const VectorOfElements & a, const VectorOfElements & b)
	{
		if (a.size() != b.size())
			errorMessage("min - two VectorOfElements have different sizes");		

		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = min(a[i],b[i]);
		return c;
	}

	VectorOfElements min(const VectorOfElements & a, const VectorOfElements & b, TypeID type)
	{
		return min(convert(type, a, false), convert(type, b, false));
	}	
	
	VectorOfElements minAbs(const VectorOfElements & a, const VectorOfElements & b)
	{
		if (a.size() != b.size())
			errorMessage("minAbs: two VectorOfElements have different sizes");		

		acl::TypeID type(getElementType(a));
		
		VectorOfElements c(a.size());
		copy(select(b, a, fabs(a)<fabs(b), type), c);
		return c;
	}
	
	VectorOfElements max(const VectorOfElements & a, const VectorOfElements & b)
	{
		if (a.size() != b.size())
			errorMessage("max - two VectorOfElements have different sizes");		

		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = max(a[i],b[i]);
		return c;
	}

	VectorOfElements max(const VectorOfElements & a, const VectorOfElements & b, TypeID type)
	{
		return max(convert(type, a, false), convert(type, b, false));
	}

	
	VectorOfElements copysign(const VectorOfElements & a, const VectorOfElements & b)
	{
		if (a.size() != b.size())
			errorMessage("max - two VectorOfElements have different sizes");		

		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = copysign(a[i],b[i]);
		return c;
	}

	VectorOfElements copysign(const VectorOfElements & a,
	                          const VectorOfElements & b,
	                          TypeID t)
	{
		return copysign(convert(t,a,false),convert(t,b,false));
	}

	VectorOfElements sign(const VectorOfElements & a)
	{
		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = sign(a[i]);
		return c;
	}
	
	VectorOfElements excerpt(const VectorOfElements & source, const VectorOfElements & filter)
	{
		if (filter.size() != 1)
			errorMessage("exerpt - filter has more than 1 component");		

		VectorOfElements c(source.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
		{
			c[i] = excerpt(source[i], filter[0]);
		}
		return c;
	}
		
	VectorOfElements select(const VectorOfElements & a,
	                        const VectorOfElements & b,
	                        const VectorOfElements & c)
	{ 
		VectorOfElements r(a.size());
		if (c.size() == 1)
		{
			if (a.size() == b.size())
			{
				using namespace elementOperators;
				for (unsigned int i(0); i < a.size(); ++i)
					r[i] = select(a[i], b[i], c[0]);
				return r;
			}		
			errorMessage("select - VectorOfElements a and b have different sizes");		
		}
		else
		{
			if (a.size() == b.size() && c.size() == b.size())
			{
				using namespace elementOperators;
				for (unsigned int i(0); i < c.size(); ++i)
					r[i] = select(a[i], b[i], c[i]);
				return r;
			}		
			errorMessage("select - VectorOfElements a, b and c have different sizes");		
		}
		return r;
	}

	VectorOfElements select(const VectorOfElements & a,
	                        const VectorOfElements & b,
	                        const VectorOfElements & c, 
	                        TypeID t)
	{
		return select(convert(t,a,false), 
		              convert(t,b,false), 
		              convert(TYPE_SELECT[t],c,false));
	}

	VectorOfElements select(const VectorOfElements & b,
	                        const VectorOfElements & c, 
	                        TypeID t)
	{
		return select(convert(t,generateVEConstantN(b.size(),0),false), 
		              convert(t,b,false), 
		              convert(TYPE_SELECT[t],c,false));
	}

	VectorOfElements mad(const VectorOfElements & a, 
	                     const VectorOfElements & b, 
	                     const VectorOfElements & c)
	{ 
		if (a.size() != 1 || b.size() != 1 || c.size() != 1)
			errorMessage("mad - the function defined only fo sigle value VectorOfElements");		
		using namespace elementOperators;
		VectorOfElements r(1);
		r[0] = mad(a[0], b[0], c[0]);
		return r;
	}

	VectorOfElements mad(const VectorOfElements & a, 
	                     const VectorOfElements & b, 
	                     const VectorOfElements & c,
	                     TypeID t)
	{ 
		return mad(convert(t,a,false), 
		              convert(t,b,false), 
		              convert(t,c,false));
	}
	
		
	VectorOfElements powI(const VectorOfElements & a, unsigned int i)
	{ 
		VectorOfElements c(1);
		if (a.size() != 1)
			errorMessage("powI - the input vector has number of elements unequal to 1");		

		using namespace elementOperators;
		c[0] = powI(a[0], i);
		return c;
	}

	VectorOfElements exp(const VectorOfElements & a)
	{ 
		if (a.size() != 1)
			errorMessage("exp- the input vector has number of elements unequal to 1");		

		VectorOfElements c(1);
		c[0] = elementOperators::exp(a[0]);
		return c;
	}

	VectorOfElements sqrt(const VectorOfElements & a)
	{ 
		if (a.size() != 1)
			errorMessage("sqrt- the input vector has number of elements unequal to 1");		

		VectorOfElements c(1);
		c[0] = elementOperators::sqrt(a[0]);
		return c;
	}

	VectorOfElements rsqrt(const VectorOfElements & a)
	{ 
		if (a.size() != 1)
			errorMessage("rsqrt- the input vector has number of elements unequal to 1");		

		VectorOfElements c(1);
		c[0] = elementOperators::rsqrt(a[0]);
		return c;
	}
	
	VectorOfElements fabs(const VectorOfElements & a)
	{ 
		VectorOfElements c(a.size());
		for(unsigned int i(0); i < a.size(); ++i)
			c[i] = elementOperators::fabs(a[i]);
		return c;
	}

	VectorOfElements abs(const VectorOfElements & a)
	{ 
		VectorOfElements c(a.size());
		for(unsigned int i(0); i < a.size(); ++i)
			c[i] = elementOperators::abs(a[i]);
		return c;
	}

	VectorOfElements abs_diff(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		if (a.size() != b.size())
			errorMessage("abs_diff - two VectorOfElements have different sizes");		
		
		VectorOfElements c(a.size());
		for(unsigned int i(0); i < a.size(); ++i)
			c[i] = elementOperators::abs_diff(a[i],b[i]);
		return c;
	}
	
	VectorOfElements floor(const VectorOfElements & a)
	{ 
		VectorOfElements c(a.size());
		for(unsigned int i(0); i < a.size(); ++i)
			c[i] = elementOperators::floor(a[i]);
		return c;
	}
	
	VectorOfElements convert(acl::TypeID type,const VectorOfElements & a, bool strong)
	{
		VectorOfElements c(a.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < c.size(); ++i)
			c[i] = convert(type,a[i],strong);
		return c;
	}
	
	
		
		
	VectorOfElements log(const VectorOfElements & a)
	{ 
		if (a.size() != 1)
			errorMessage("log - the input vector has number of elements unequal to 1");		
		VectorOfElements c(1);
		c[0] = elementOperators::log(a[0]);		
		return c;
	}


	VectorOfElements log10(const VectorOfElements & a)
	{ 
		if (a.size() != 1)
			errorMessage("log10 - the input vector has number of elements unequal to 1");		
		VectorOfElements c(1);
		c[0] = elementOperators::log10(a[0]);		
		return c;
	}
		

	VectorOfElements minElement(const VectorOfElements & a)
	{
		VectorOfElements c(1);
		c[0]=min(a[0], a[1]);
		for(unsigned int i(2); i<a.size(); ++i)
			c[0]=acl::elementOperators::min(c[0], a[i]);
		return c;
	}

	VectorOfElements minAbsElement(const VectorOfElements & a)
	{
		VectorOfElements c(1);
		acl::TypeID type(getElementType(a));
		using namespace elementOperators;
		c[0]=select(a[0], a[1], convert(TYPE_SELECT[type],
		                                fabs(a[1])<fabs(a[0]),
		                                false));
		for(unsigned int i(2); i<a.size(); ++i)
			c[0]=select(c[0], a[i], convert(TYPE_SELECT[type],
			                                fabs(a[i])<fabs(c[0]), 
			                                false));
		return c;
	}
	
	VectorOfElements maxElement(const VectorOfElements & a)
	{
		VectorOfElements c(1);
		c[0]=min(a[0], a[1]);
		for(unsigned int i(2); i < a.size(); ++i)
			c[0]=acl::elementOperators::max(c[0], a[i]);
		return c;
	}

	VectorOfElements sumOfElements(const VectorOfElements & a)
	{
		VectorOfElements c(1);
		c[0]=acl::elementOperators::operator+(a[0], a[1]);
		for(unsigned int i(2); i < a.size(); ++i)
			c[0]=acl::elementOperators::operator+(c[0], a[i]);
		return c;
	}

	VectorOfElements andOfElements(const VectorOfElements & a)
	{
		VectorOfElements c(1);
		c[0]=acl::elementOperators::operator&&(a[0], a[1]);
		for(unsigned int i(2); i < a.size(); ++i)
			c[0]=acl::elementOperators::operator&&(c[0], a[i]);
		return c;
	}

	VectorOfElements orOfElements(const VectorOfElements & a)
	{
		VectorOfElements c(1);
		c[0]=acl::elementOperators::operator||(a[0], a[1]);
		for(unsigned int i(2); i < a.size(); ++i)
			c[0]=acl::elementOperators::operator||(c[0], a[i]);
		return c;
	}
	
	VectorOfElements cat(const VectorOfElements & a, const VectorOfElements & b)
	{ 
		VectorOfElements c(a.size() + b.size());
		using namespace elementOperators;			
		for (unsigned int i(0); i < a.size(); ++i)
			c[i] = ( a[i]);
		for (unsigned int i(0); i < b.size(); ++i)
			c[a.size() + i] = b[i];
		return c;
	}

	VectorOfElementsData cat(const VectorOfElementsData & a, const VectorOfElementsData & b)
	{ 
		VectorOfElementsData c(a.size() + b.size());
		using namespace elementOperators;			
		for (unsigned int i(0); i < a.size(); ++i)
			c[i] = ( a[i]);
		for (unsigned int i(0); i < b.size(); ++i)
			c[a.size() + i] = b[i];
		return c;
	}

	VectorOfElements cat(const VectorOfElements & a, 
	                     const VectorOfElements & b, 
	                     const VectorOfElements & c)
	{ 
		VectorOfElements z(cat(a,b));
		return cat(z,c);
	}

	template <typename T> VectorOfElements operator+=(const VectorOfElements & a, const T & b)
	{
		return a+=generateVEConstant(b);
	}

	template VectorOfElements operator+=(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator+=(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator+=(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator+=(const VectorOfElements & a, const cl_int & b);	

	template <typename T> VectorOfElements operator-=(const VectorOfElements & a, const T & b)
	{
		return a-=generateVEConstant(b);
	}

	template VectorOfElements operator-=(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator-=(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator-=(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator-=(const VectorOfElements & a, const cl_int & b);	

	template <typename T> VectorOfElements operator*=(const VectorOfElements & a, const T & b)
	{
		return a*=generateVEConstant(b);
	}

	template VectorOfElements operator*=(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator*=(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator*=(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator*=(const VectorOfElements & a, const cl_int & b);	

	template <typename T> VectorOfElements operator/=(const VectorOfElements & a, const T & b)
	{
		return a/=generateVEConstant(b);
	}

	template VectorOfElements operator/=(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator/=(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator/=(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator/=(const VectorOfElements & a, const cl_int & b);	

	template <typename T> VectorOfElements operator+(const VectorOfElements & a, const T & b)
	{
		return a+generateVEConstant(b);
	}

	template VectorOfElements operator+(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator+(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator+(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator+(const VectorOfElements & a, const cl_int & b);	
	template VectorOfElements operator+(const VectorOfElements & a, const AVec<cl_double> & b);
	template VectorOfElements operator+(const VectorOfElements & a, const AVec<cl_float> & b);
	template VectorOfElements operator+(const VectorOfElements & a, const AVec<cl_uint> & b);
	template VectorOfElements operator+(const VectorOfElements & a, const AVec<cl_int> & b);	
	
	template <typename T> inline VectorOfElements operator+(const T & a, const VectorOfElements & b)
	{
		return generateVEConstant(a)+b;
	}

	template VectorOfElements operator+(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator+(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator+(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator+(const cl_uint & a, const VectorOfElements & b);	
	template VectorOfElements operator+(const AVec<cl_double> & a, const VectorOfElements & b);
	template VectorOfElements operator+(const AVec<cl_float> & a, const VectorOfElements & b);
	template VectorOfElements operator+(const AVec<cl_int> & a, const VectorOfElements & b);
	template VectorOfElements operator+(const AVec<cl_uint> & a, const VectorOfElements & b);	
	
	template <typename T> VectorOfElements operator-(const VectorOfElements & a, const T & b)
	{
		return a-generateVEConstant(b);
	}

	template VectorOfElements operator-(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator-(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator-(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator-(const VectorOfElements & a, const cl_int & b);	
	template VectorOfElements operator-(const VectorOfElements & a, const AVec<cl_double> & b);
	template VectorOfElements operator-(const VectorOfElements & a, const AVec<cl_float> & b);
	template VectorOfElements operator-(const VectorOfElements & a, const AVec<cl_uint> & b);
	template VectorOfElements operator-(const VectorOfElements & a, const AVec<cl_int> & b);	
	
	template <typename T> VectorOfElements operator-(const T & a, const VectorOfElements & b)
	{
		return generateVEConstant(a)-b;
	}

	template VectorOfElements operator-(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator-(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator-(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator-(const cl_uint & a, const VectorOfElements & b);	
	template VectorOfElements operator-(const AVec<cl_double> & a, const VectorOfElements & b);
	template VectorOfElements operator-(const AVec<cl_float> & a, const VectorOfElements & b);
	template VectorOfElements operator-(const AVec<cl_int> & a, const VectorOfElements & b);
	template VectorOfElements operator-(const AVec<cl_uint> & a, const VectorOfElements & b);	

	template <typename T> VectorOfElements operator*(const VectorOfElements & a,const T & b)
	{
		return a*generateVEConstant(b);
	}

	template VectorOfElements operator*(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator*(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator*(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator*(const VectorOfElements & a, const cl_int & b);	
	template VectorOfElements operator*(const VectorOfElements & a, const AVec<cl_double> & b);		
	template VectorOfElements operator*(const VectorOfElements & a, const AVec<cl_float> & b);		
	template VectorOfElements operator*(const VectorOfElements & a, const AVec<cl_int> & b);		

	template <typename T> VectorOfElements operator*(const T & a, const VectorOfElements & b)
	{
		return generateVEConstant(a)*b;
	}

	template VectorOfElements operator*(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator*(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator*(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator*(const cl_uint & a, const VectorOfElements & b);	

	template <typename T> VectorOfElements operator/(const VectorOfElements & a, const T & b)
	{
		return a/generateVEConstant(b);
	}

	template VectorOfElements operator/(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator/(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator/(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator/(const VectorOfElements & a, const cl_int & b);	

	template <typename T> VectorOfElements operator/(const T & a, const VectorOfElements & b)
	{
		return generateVEConstant(a)/b;
	}

	template VectorOfElements operator/(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator/(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator/(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator/(const cl_uint & a, const VectorOfElements & b);	

	template <typename T> VectorOfElements operator%(const VectorOfElements & a, const T & b)
	{
		return a%generateVEConstant(b);
	}

	template VectorOfElements operator%(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator%(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator%(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator%(const VectorOfElements & a, const cl_int & b);	

	template <typename T> VectorOfElements operator%(const T & a, const VectorOfElements & b)
	{
		return generateVEConstant(a)%b;
	}

	template VectorOfElements operator%(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator%(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator%(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator%(const cl_uint & a, const VectorOfElements & b);	
	
	template <typename T> VectorOfElements operator>(const VectorOfElements & a, const T & b)
	{
		return a>generateVEConstant(b);
	}

	template VectorOfElements operator>(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator>(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator>(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator>(const VectorOfElements & a, const cl_int & b);	
	
	template <typename T> VectorOfElements operator>(const T & a, const VectorOfElements & b)
	{
		return generateVEConstant(a)>b;
	}

	template VectorOfElements operator>(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator>(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator>(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator>(const cl_uint & a, const VectorOfElements & b);	

	template <typename T> VectorOfElements operator<(const VectorOfElements & a, const T & b)
	{
		return a<generateVEConstant(b);
	}

	template VectorOfElements operator<(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator<(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator<(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator<(const VectorOfElements & a, const cl_int & b);	
	
	template <typename T> VectorOfElements operator<(const T & a, const VectorOfElements & b)
	{
		return generateVEConstant(a)<b;
	}

	template VectorOfElements operator<(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator<(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator<(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator<(const cl_uint & a, const VectorOfElements & b);	

	template <typename T> VectorOfElements operator>=(const VectorOfElements & a, const T & b)
	{
		return a>=generateVEConstant(b);
	}

	template VectorOfElements operator>=(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator>=(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator>=(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator>=(const VectorOfElements & a, const cl_int & b);	
	
	template <typename T> VectorOfElements operator>=(const T & a, const VectorOfElements & b)
	{
		return generateVEConstant(a)>=b;
	}

	template VectorOfElements operator>=(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator>=(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator>=(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator>=(const cl_uint & a, const VectorOfElements & b);	

	template <typename T> VectorOfElements operator<=(const VectorOfElements & a, const T & b)
	{
		return a<=generateVEConstant(b);
	}

	template VectorOfElements operator<=(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator<=(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator<=(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator<=(const VectorOfElements & a, const cl_int & b);	
	
	template <typename T> VectorOfElements operator<=(const T & a, const VectorOfElements & b)
	{
		return generateVEConstant(a)<=b;
	}

	template VectorOfElements operator<=(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator<=(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator<=(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator<=(const cl_uint & a, const VectorOfElements & b);	

	template <typename T> VectorOfElements operator==(const VectorOfElements & a,
	                                                  const T & b)
	{
		return a == generateVEConstant(b);
	}

	template VectorOfElements operator==(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator==(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator==(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator==(const VectorOfElements & a, const cl_int & b);	
	

	template <typename T> VectorOfElements operator==(const T & a,
	                                                  const VectorOfElements & b)
	{
		return generateVEConstant(a) == b;
	}

	template VectorOfElements operator==(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator==(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator==(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator==(const cl_uint & a, const VectorOfElements & b);	

	template <typename T> VectorOfElements operator!=(const VectorOfElements & a, const T & b)
	{
		return a!=generateVEConstant(b);
	}

	template VectorOfElements operator!=(const VectorOfElements & a, const cl_double & b);
	template VectorOfElements operator!=(const VectorOfElements & a, const cl_float & b);
	template VectorOfElements operator!=(const VectorOfElements & a, const cl_uint & b);
	template VectorOfElements operator!=(const VectorOfElements & a, const cl_int & b);	
	
	template <typename T> VectorOfElements operator!=(const T & a, const VectorOfElements & b)
	{
		return generateVEConstant(a)!=b;
	}

	template VectorOfElements operator!=(const cl_double & a, const VectorOfElements & b);
	template VectorOfElements operator!=(const cl_float & a, const VectorOfElements & b);
	template VectorOfElements operator!=(const cl_int & a, const VectorOfElements & b);
	template VectorOfElements operator!=(const cl_uint & a, const VectorOfElements & b);	

	vector<Element> gcLength2(const VectorOfElements & a, const VectorOfElements & l2)
	{
		auto type(getElementType(l2));
		auto temp(acl::generateVEPrivateVariable(1u, type));
		vector<Element> code;
		code<<(temp=subVE(a,0));
		code<<(l2=temp*temp);
		for(unsigned int i(1); i<a.size(); ++i)
		{
			code<<(temp=subVE(a,i));
			code<<(l2+=temp*temp);
		}
		return code;
	}

	vector<Element> gcLength(const VectorOfElements & a, const VectorOfElements & l)
	{
		auto code(gcLength2(a,l));
		code.push_back((l=sqrt(l))[0]);
		return code;
	}
	
	vector<Element> gcNormalize(const VectorOfElements & a)
	{
		auto type(getElementType(a));
		auto len(acl::generateVEPrivateVariable(1u, type));
		auto code(gcLength(a,len));
		auto res(a /= len);
		copy(res.begin(), res.end(), back_inserter(code));
			
		return code;
	}
	

	
} // namespace acl
