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


/// \file aslVectors.h definition of class АVec<T>

#ifndef ASLVECTORSDYNAMICLENGTHOPERATIONS_H
#define ASLVECTORSDYNAMICLENGTHOPERATIONS_H


#include "../aslUtilities.h"
#include <cmath>


namespace asl
{
	
	/// The function checks whether the sizes are equal \relates AVec
	template <typename T1, typename T2> 
		inline bool compatibleSizes(AVec<T1> a, AVec<T2> b);

	/// \relates AVec	
	template <typename T> 
		inline const T l2(const AVec<T> & a);

	/// \relates AVec	
	template <typename T> 
		inline const AVec<> normalize(const AVec<T> & a);
	
	/// \relates AVec	
	template <typename T> 
		inline const AVec<T> operator+(const AVec<T> & a, const AVec<T> & b);
	/// \relates AVec
	template <typename T> 
		inline const AVec<T> operator-(const AVec<T> & a, const AVec<T> & b);
	/// \relates AVec
	template <typename T> 
		inline const AVec<T> operator-(const AVec<T> & a);
	/// \relates AVec
	template <typename T> 
		inline const T operator*(const AVec<T> & a, const AVec<T> & b);	
	/// \relates AVec
	template <typename T> 
		inline const AVec<T> operator*(const T & a, const AVec<T> & b);	
	/// \relates AVec
	template <typename T> 
		inline const AVec<T> operator*(const AVec<T> & a, const T & b);		
	/// \relates AVec
	template <typename T> 
		inline const AVec<T> operator/(const AVec<T> & a, const T & b);		

	/// \relates AVec
	template <typename T> inline const AVec<T> & operator+=(AVec<T> & a,
	                                                        const AVec<T> & b);
	/// \relates AVec
	template <typename T> inline const AVec<T> & operator-=(AVec<T> & a,
	                                                        const AVec<T> & b);
	/// \relates AVec
	template <typename T> inline const AVec<T> & operator*=(AVec<T> & a,
	                                                        const T & b);		

	/// \relates AVec
	template <typename T> inline const bool operator==(const AVec<T> & a,
	                                                   const AVec<T> & b);
	

	/// \relates AVec
	template <typename T> inline const bool operator!=(const AVec<T> & a,
	                                                   const AVec<T> & b);
	

	/// \relates AVec
	template <typename T> 
		inline const AVec<T> crossProduct(const AVec<T> & a, const AVec<T> & b);	

	/// returns minimal component \relates AVec 
	template <typename T> inline const T minComponent(const AVec<T> & a);		
	/// returns maximal component  \relates AVec
	template <typename T> inline const T maxComponent(const AVec<T> & a);		

	/// returns summ of all components \relates AVec
	template <typename T> inline T sumOfElements(const AVec<T> & a);
	/// returns product of all components \relates AVec
	template <typename T> inline T productOfElements(const AVec<T> & a);

	/// returns vector which elements are product of corresponding elements of \p a and \p b \relates AVec
	template <typename T> inline const AVec<T> productOfElements(const AVec<T> & a,
	                                                             const AVec<T> & b);
	/// returns vector which elements are division of corresponding elements of \p a and \p b \relates AVec
	template <typename T> inline const AVec<T> divisionOfElements(const AVec<T> & a,
	                                                              const AVec<T> & b);

	/// returns \p true in case when all components of \p a more or  then 0 \relates AVec 
	template <typename T> inline const bool positive(const AVec<T> & a);

	/// returns \p true in case when all components of \p a more or equal 0 \relates AVec 
	template <typename T> inline const bool nonNegative(const AVec<T> & a);
	
	/// returns \p true in case when all components of \p a more then 0 \relates AVec 
	inline const AVec<int> floor(const AVec<> & a);

	/// returns \p true in case when all components of \p a more then 0 \relates AVec 
	inline const AVec<int> round(const AVec<> & a);
	
	/// computes polynom for \p x with \p coefs 
	/**
		 \relates AVec
		 The polinom is \f$ x^{n-1}*coefs_0+x^{n-2}coefs_1+...+coefs_{n-1} \f$
	*/
	inline double computePolynom(double x, AVec<> &coefs);
	
	/// \relates AVec
	template <typename T> 
		inline std::ostream & operator<<(std::ostream & output, const AVec<T> & a);

	/// returns \p true in case when all components of \p a more then 0 \relates AVec 
	inline const AVec<> swapXZ(const AVec<> & a);
	
//---------------- Implementation ----------------

	template <typename T1, typename T2> 
	inline bool compatibleSizes(AVec<T1> a, AVec<T2> b)
	{
		return a.getSize() == b.getSize();
	}


	template <typename T> inline const T l2(const AVec<T> & a)
	{
		return a * a;
	}

	template <typename T> inline const AVec<> normalize(const AVec<T> & a)
	{
		return AVec<>(a) / sqrt(l2(a));
	}
	
	template <typename T> inline const AVec<T> operator+(const AVec<T> & a,
	                                                     const AVec<T> & b){
		if (!compatibleSizes(a, b))
			errorMessage("(AVec; operator+) Vector sizes are incompatible");
		AVec<T> c(a.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i) 
			c[i] = a[i] + b[i];
		return c;
	}

	template <typename T> inline const AVec<T> operator-(const AVec<T> & a, const AVec<T> & b)
	{
		if (!compatibleSizes (a,b))
			errorMessage("(AVec; operator-) Vector sizes are incompatible");
		AVec<T> c(a.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i) 
			c[i] =a[i]-b[i];
		return c;
	}

	template <typename T> inline const AVec<T> operator-(const AVec<T> & a)
	{
		AVec<T> c(a.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i) 
			c[i] =-a[i];
		return c;
	}
	
	template <typename T> inline const T operator*(const AVec<T> & a, const AVec<T> & b){
		if (!compatibleSizes (a,b))
			errorMessage("(AVec; operator*) Vector sizes are incompatible");
		T s(0);
		for (unsigned int i(0); i < a.getSize(); ++i) 
			s+=a[i]*b[i];
		return s;
	}

	template <typename T> inline const AVec<T> operator*(const AVec<T> & a, const T & b){
		AVec<T> c(a.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i) 
			c[i] =a[i]*b;
		return c;
	}

	template <typename T> inline const AVec<T> operator*(const T & a, const AVec<T> & b){
		return b*a;
	}

	template <typename T> inline const AVec<T> operator/(const AVec<T> & a, const T & b){
		AVec<T> c(a.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i) 
			c[i] =a[i]/b;
		return c;
	}
		
	template <typename T> inline const AVec<T> & operator+=(AVec<T> & a, const AVec<T> & b){
		if (!compatibleSizes (a,b))
			errorMessage("Vector sizes are incompatible");
		for (unsigned int i(0); i < a.getSize(); ++i) 
			a[i]+=b[i];
		return a;
	}

	template <typename T> inline const AVec<T> & operator-=(AVec<T> & a, const AVec<T> & b){
		if (!compatibleSizes (a,b))
			errorMessage("Vector sizes are incompatible");
		for (unsigned int i(0); i < a.size; ++i) 
			a[i]-=b[i];
		return a;
	}

	template <typename T> inline const AVec<T> & operator*=(AVec<T> & a, const T & b){
		for (unsigned int i(0); i < a.getSize(); ++i) 
			a[i] *= b;
		return a;
	}	

	template <typename T> inline const bool operator==(const AVec<T> & a,
	                                                   const AVec<T> & b)
	{
		if (!compatibleSizes (a, b))
			return false;
		bool c(true);
		for (unsigned int i(0); i < a.getSize(); ++i) 
			c = c && (a[i] == b[i]);
		return c;
	}


	template <typename T> inline const bool operator!=(const AVec<T> & a,
	                                                   const AVec<T> & b)
	{
		return !(a == b);
	}

	template <typename T> 
		inline const AVec<T> crossProduct(const AVec<T> & a, const AVec<T> & b)
	{
		if (!compatibleSizes (a,b))
			errorMessage("(AVec; crossProduct) Vector sizes are incompatible");
		if (a.getSize()>3)
			errorMessage("(AVec; crossProduct) number of components is more than 3");
		if (a.getSize()<2)
			errorMessage("(AVec; crossProduct) number of components is less than 2");

		AVec<T> res(1);
		if(a.getSize() == 2)
		{
			res[0] = a[0]*b[1]-a[1]*b[0];
		}
		if(a.getSize() == 3)
		{
			res.resize(3);
			res[0] = a[1]*b[2]-a[2]*b[1];
			res[1] = a[2]*b[0]-a[0]*b[2];
			res[2] = a[0]*b[1]-a[1]*b[0];
		}
		return res;
	}
	
	template <typename T> inline const T minComponent(const AVec<T> & a)
	{
		T ma(a[0]);
		for (unsigned int i(1); i < a.getSize(); ++i) 
			ma = std::min(ma, a[i]);
		return ma;
	}
		
	template <typename T> inline const T maxComponent(const AVec<T> & a)
	{
		if (a.getSize()<1) errorMessage("Vector size less than 1");
		T ma(a[0]);
		for (unsigned int i(1); i < a.getSize(); ++i) 
			ma=std::max(ma,a[i]);
		return ma;
	}

	template <typename T> inline AVec<T> subAVec(const AVec<T> & source,
	                       unsigned int start, unsigned int end)
	{
		if (source.getSize() <= end )
			errorMessage("subAVec: attempt to copy besides the vector range");
		
		AVec<T> destination(1 + end - start);
		for (unsigned int i(start); i <= end; ++i)
			destination[i - start] = source[i];
		return destination;
	}
	
	template <typename T> inline T sumOfElements(const AVec<T> & a)
	{
		T s(0);
		for (unsigned int i(0); i < a.getSize(); ++i)
			s += a[i];
		return s;
	}
		
	template <typename T> inline T productOfElements(const AVec<T> & a)
	{
		T p(1);
		for (unsigned int i(0); i < a.getSize(); ++i)
			p *= a[i];
		return p;
	}

	template <typename T> inline const AVec<T> productOfElements(const AVec<T> & a,
	                                                             const AVec<T> & b)
	{
		if (!compatibleSizes (a, b))
			errorMessage("(AVec; productOfElements) Vector sizes are incompatible");
		AVec<T> c(a.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i)
			c[i] = a[i] * b[i];
		return c;
	}

	template <typename T> inline const AVec<T> divisionOfElements(const AVec<T> & a, const AVec<T> & b)
	{
		if (!compatibleSizes (a, b))
			errorMessage("(AVec; divisionOfElements) Vector sizes are incompatible");

		AVec<T> c(a.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i)
			c[i] = a[i] / b[i];

		return c;
	}

	template <typename T> inline const bool positive(const AVec<T> & a)
	{
		if (!a.getSize())
			errorMessage("(AVec; positive) Vector size is zero");
		bool res(a[0]>0);
		for (unsigned int i(1); i < a.getSize(); ++i)
			res &= a[i]>0;
		return res;	
	}

	template <typename T> inline const bool nonNegative(const AVec<T> & a)
	{
		if (!a.getSize())
			errorMessage("(AVec; positive) Vector size is zero");
		bool res(a[0]>=0);
		for (unsigned int i(1); i < a.getSize(); ++i)
			res &= a[i]>=0;
		return res;	
	}
	
	inline const AVec<int> floor(const AVec<> & a)
	{
		if (!a.getSize())
			errorMessage("(AVec; floor) Vector size is zero");
		AVec<int> res(a.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i)
			res[i] = std::floor(a[i]);
		return res;	
	}

	inline const AVec<int> round(const AVec<> & a)
	{
		if (!a.getSize())
			errorMessage("(AVec; round) Vector size is zero");
		AVec<int> res(a.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i)
			res[i] = std::round(a[i]);
		return res;	
	}

	
	inline double computePolynom (double x, AVec<> &coefs)
	{	
		if (coefs.getSize() < 1)
			errorMessage("Error: (asl::computePolynom) size of \"coefs\" less than 1");
		double p;
		p=coefs[0];
		for (unsigned int i(1); i < coefs.getSize(); ++i)
			p=p*x+coefs[i];
		return p;	
	}

	/// returns \p true in case when all components of \p a more then 0 \relates AVec 
	inline const AVec<> swapXZ(const AVec<> & a)
	{
		if (a.getSize()<3)
			errorMessage("(AVec; swapXZ) Vector size less than 3");
		AVec<> res(a);
		std::swap(res[0],res[2]);
		return res;		
	}
	
} // asl

#endif //ASLVECTORSDYNAMICLENGTH_H
