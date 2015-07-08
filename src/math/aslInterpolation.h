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

#ifndef ASLINTERPOLATION
#define ASLINTERPOLATION


#include "../aslUtilities.h"
#include <math.h>


namespace asl
{
	class ElementBase;
	typedef std::shared_ptr<ElementBase> Element;

	template <class Func> class UniversalFunction
	{
		public:
			inline double operator()(double x);
			inline double operator()(Element x);
	};

	///Linear spline function
	/**
	 \ingroup Splines

	this class defines a spline function \f$ (1-|x/r_0|),~~~ |x|<r_0 \f$
	 \image html linearS.png ""
	*/
	class LinearSpline
	{
		public:
			template <typename T> inline T operator(T x,T r0)
			{
				return fabs(x)<r0 ? 1.- fabs(x)/r0 : 0;
			}
	};

	///Quadratic spline function
	/**
	 \ingroup Splines

	this class defines a spline function \f$ (1-|x|/r_0)^2,~~~ |x|<r_0 \f$
	 \image html quadraticS.png ""
	*/
	class QuadraticSpline
	{
		public:
			template <typename T> inline T operator(T x,T r0)
			{
				T a(fabs(x));
				T b((1.- a/r0));
				return a<r0 ? b*b : 0;
			}
	};

	///Qubic spline function
	/**	
	 \ingroup Splines

	this class defines a spline function \f$ (1-|x|/r_0)^3,~~~ |x|<r_0 \f$	 
	 \image html qubicS.png ""
	*/
	class QubicSpline
	{
		public:
			template <typename T> inline T operator(T x,T r0)
			{
				T a(fabs(x));
				T b((1.- a/r0));
				return a<r0 ? b*b*b : 0;
			}
	};

	///Qubic spline function
	/**	
	 \ingroup Splines

	this class defines a spline function \f$ (1-|x|/r_0)^3,~~~ |x|<r_0 \f$	 
	 \image html qubicS1.png ""
	*/
	class QubicSpline1
	{
		public:
			template <typename T> inline T operator(T x,T r0)
			{
				T a(fabs(x));
				T b(a/r0);
				T b2(b*b);
				return b < 1 ? 2. * b2*b - 3. * b2 + 1 : 0;
			}
	};
	
// --------------------------- Implementation ---------------------

	template <class Func> inline double UniversalFunction<Func>::operator()(double x)
	{
		return Func(x);
	}

	template <class Func>  inline Element UniversalFunction<Func>::operator()(Element x)
	{
		return Func(x);
	}
		
} // asl

#endif

