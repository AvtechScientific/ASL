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

#ifndef ASLVECTORSDYNAMICLENGTH_H
#define ASLVECTORSDYNAMICLENGTH_H


#include "../aslUtilities.h"
#include <math.h>


namespace asl
{
	
	/// class algebraic vector. 
	/// The class is an implementation of a dynamic array with defined algebraic operations
	template <typename T = double> class AVec
	{
		private:
			T* x;	
			unsigned int size;

		public:
			typedef T Type;
			inline AVec();						
			inline explicit AVec(unsigned int s);			
			inline AVec(unsigned int s, T a);			
			inline ~AVec();
			template <typename Tv> 
				inline explicit AVec(const AVec<Tv> & a);
			inline AVec(const AVec<T> & a);
			template <typename Tv> 
				inline const AVec<T> & operator=(const AVec<Tv> &a);
			template <typename Tv> 
				inline const AVec<T> & operator=(const std::vector<Tv> &a);
			inline const AVec<T> & operator=(const AVec<T> &a);
			inline T& operator[](unsigned int i);
			inline const T& operator[](unsigned int i) const;
			inline const unsigned int & getSize() const;
			inline void resize(unsigned int newSize);
	};

	/// \relates AVec
	template <typename T> 
		inline const unsigned int nD(const AVec<T> a);
	
	/// \relates AVec
	template <typename T> inline AVec<T> makeAVec(T a1);
	/// \relates AVec
	template <typename T> inline AVec<T> makeAVec(T a1, T a2);
	/// \relates AVec
	template <typename T> inline AVec<T> makeAVec(T a1, T a2, T a3);
	/// \relates AVec
	template <typename T> inline AVec<T> makeAVec(T a1, T a2, T a3, T a4);
	
	/// \relates AVec
	template <typename T> 
		inline std::ostream & operator<<(std::ostream & output, const AVec<T> & a);

	
//---------------- Implementation ----------------
	template <typename T> inline AVec<T>::AVec():
		x(new T[1]),
		size(1)
	{
		x[0] =0;
	}

	template <typename T> inline AVec<T>::AVec(unsigned int s):
		x(new T[s]),
		size(s)
	{
		if (size<1) errorMessage("Vector size less than 1");
		memset(x, 0, sizeof(T)*s);
	}

	template <typename T> inline AVec<T>::AVec(unsigned int s, T a):
		x(new T[s]),
		size(s)
	{
		if (size<1) errorMessage("Vector size less than 1");		
		for (unsigned int i(0); i<size; ++i)
			x[i] =a;
	}
		
	template <typename T> inline AVec<T>::~AVec()
	{	
		delete[] x;
	}

	template <typename T> template<typename Tv> inline AVec<T>::AVec(const AVec<Tv> & a):
		x(new T[a.getSize()]),
		size(a.getSize())	
	{
		for (unsigned int i(0); i < size; ++i)
			x[i] = a[i];		
	}

	template <typename T> inline AVec<T>::AVec(const AVec<T> & a):
		x(new T[a.getSize()]),
		size(a.getSize())	
	{
		for (unsigned int i(0); i < size; ++i)
			x[i] = a[i];		
	}

		
	template <typename T> template<typename Tv> 
		inline const AVec<T> & AVec<T>::operator=(const AVec<Tv> & a)
	{
		resize(a.getSize());
		for (unsigned int i(0); i<size; ++i)
			x[i] =a[i];	
		return *this;
	}
		
	template <typename T>
		inline const AVec<T> & AVec<T>::operator=(const AVec<T> & a)
	{
		resize(a.getSize());
		for (unsigned int i(0); i<size; ++i)
			x[i] =a[i];		
		return *this;
	}

	template <typename T> template<typename Tv> 
		inline const AVec<T> & AVec<T>::operator=(const std::vector<Tv> & a)
	{
		resize(a.size());
		for (unsigned int i(0); i<size; ++i)
			x[i] =a[i];	
		return *this;
	}

		
	template <typename T> 
	inline T& AVec<T>::operator[](unsigned int i)
	{
		return x[i];
	}

	template <typename T> 
	inline const T& AVec<T>::operator[](unsigned int i) const
	{
		return x[i];
	}

	template <typename T> 
	inline const unsigned int & AVec<T>::getSize() const
	{
		return size;
	}
								
	template <typename T> inline AVec<T> makeAVec(T a1)
	{
		return AVec<T>(1, a1);
	}
		
	template <typename T> inline AVec<T> makeAVec(T a1, T a2)
	{
		AVec<T> v(2); 
		v[0] = a1;
		v[1] = a2;
		return v;
	}

	template <typename T> inline AVec<T> makeAVec(T a1, T a2, T a3)
	{
		AVec<T> v(3); 
		v[0] = a1;
		v[1] = a2;
		v[2] = a3;
		return v;
	}

	template <typename T> inline AVec<T> makeAVec(T a1, T a2, T a3, T a4)
	{
		AVec<T> v(4); 
		v[0] = a1;
		v[1] = a2;
		v[2] = a3;
		v[3] = a4;
		return v;
	}
		
	template <typename T> inline std::ostream & operator<<(std::ostream & output, 
	                                                       const AVec<T> & a)
	{
		for (unsigned int i(0); i<a.getSize(); ++i)
			output << a[i]<< " ";
		return output;
	}

	template <typename T> 
		inline const unsigned int nD(const AVec<T> a)
	{
		return a.getSize();
	}
		
	template <typename T>
		inline void AVec<T>::resize(unsigned int newSize)
	{
		if (size != newSize){
			delete[] x;
			x=new T[newSize];
			size=newSize;
		}
	}
		
} // asl

#endif //ASLVECTORSDYNAMICLENGTH_H

