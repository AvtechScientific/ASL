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


/// \file aslUtilities.h useful common utilities

#ifndef ASLUTILITIES_H
#define ASLUTILITIES_H

#include <acl/aclStdIncludes.h>
#include <CL/cl.hpp>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <list>
#include <typeinfo>

/// Advanced Simulation Library
namespace asl
{
	
	/// Converts numbers or another type to string \ingroup Utilities
	template <typename T> inline std::string numToStr(T i)
	{
		std::stringstream s;
		s << i;
		return s.str();
	}


	/// Converts numbers or another type to string with given value of positions with 0 before \ingroup Utilities
	template <typename T> std::string numToStr(T i, int numberOf0)
	{
		std::stringstream s;
		s << std::setfill('0') << std::setw(numberOf0) << i;
		return s.str();
	}


	/// Converts string to number, exits if not able to convert
	template <typename T> T strToNum(std::string s);
	
	/// This class is used in order to overload function according to an integer parameter
	template <int I>class I2T{public: static const int i = I;};

	/// Realization of \f$ a^I \f$ 
	template <int I> inline const double P(const double& a){return P<I-1>(a)*a;}
	template <int I> inline const float P(const float& a){return P<I-1>(a)*a;}
	template <int I> inline const int P(const int& a){return P<I-1>(a)*a;}	

	/// \cond NEVER
	template <> inline const double P<1>(const double& a){return a;}
	template <> inline const float P<1>(const float& a){return a;}
	template <> inline const int P<1>(const int& a){return a;}	
	/// \endcond


	/// Checks the belonging to a closed interval [x1,x2], \f$ xx \in [x1,x2] \f$ 
	template <typename T> inline  bool in(const T& xx, const T& x1, const T& x2)
	{
		return (((x1-xx)*(x2-xx))<=0);
	}

	/// Checks the belonging to an open interval (x1,x2), \f$ xx \in (x1,x2) \f$
	template <typename T> inline  bool inO(const T& xx, const T& x1, const T& x2)
	{
		return (((x1 - xx) * (x2 - xx)) < 0);
	}

	static const double pi=4.*atan(1.); ///<\f$\pi\f$	

	/// Translation radians to degrees
	inline const double deg(double a){return a/pi*180.;}
	inline const float deg(float a){return a/pi*180.;}	
	/// Translation degrees to radians
	inline const double rad(double a){return a*pi/180.;}
	inline const float rad(float a){return a*pi/180.;}	

	/// Approximately equal; the precision is defined as \p p_
	inline const bool approxEqual(const double &a, const double &b, const double p_ = 1e-6)
	{
		return fabs(a-b) < p_ * fabs(a) && fabs(a-b) < p_ * fabs(b);
	}
	
	inline const bool approxEqual(const float &a, const float &b, const float p_ = 1e-6)
	{
		return fabs(a - b) < p_ * fabs(a) && fabs(a - b) < p_ * fabs(b);
	}


	/// Prints errorMessage and exits depending on the status \ingroup ErrorMessaging
	void errorMessage(cl_int status, const char *errorMessage);

	/// Prints errorMessage and exits depending on the status \ingroup ErrorMessaging
	void errorMessage(cl_int status, const std::string  & errorMessage);
	
	/// Prints " Ok" or " ERROR" depending on status \ingroup ErrorMessaging
	void errorMessage(bool status);

	/// Prints errorMessage and exits \ingroup ErrorMessaging
	void errorMessage(const char *errorMessage);

	/// Prints errorMessage and exits \ingroup ErrorMessaging
	void errorMessage(const std::string & errorMessage);


	/// Prints warningMessage \ingroup ErrorMessaging
	void warningMessage(const char *warningMessage);
	

	/// Prints warningMessage \ingroup ErrorMessaging
	void warningMessage(const std::string & warningMessage);
	

	/// Returns warningMessage \ingroup ErrorMessaging
	std::string warningString(const char *warningMessage);


	/// Prints elements of the vector separated by space
	template <typename T> std::ostream & operator<<(std::ostream & output, 
	                                                const std::vector<T> & vector);


	/// Compares two vectors
	template <typename T> bool operator==(const std::vector<T> & vector1, 
	                                      const std::vector<T> & vector2);

	
	/// Makes output 1  of \p n times 
	template <class T, int N> inline void output1OfN(const std::string &s);


	template <class T>inline void setupAll(std::vector<std::shared_ptr<T> > & v);


	template <class T>inline void computeAll(std::vector<std::shared_ptr<T> > & v);

	/// sorts two vectors with respect of the fist one
	template <class T1, class T2> inline void sortTwoVectors(std::vector<T1> & v1, std::vector<T2> & v2);

	/// reorders vector according to indeces
	template <class T> inline void reorderVector(std::vector<unsigned int> & ind, std::vector<T> & v);

	
//---------------------------- Implementations ----------------------------

	/// Prints elements of the vector separated by space
	template <typename T> std::ostream & operator<<(std::ostream & output, const std::vector<T> & vector)
	{
		typename std::vector<T>::const_iterator iter;

		for (iter = vector.begin(); iter != vector.end(); ++iter)
			output << *iter << " ";
		output << std::endl;
		return output;
	}


	/// compares two vectors
	template <typename T> bool operator==(const std::vector<T> & vector1, 
	                                      const std::vector<T> & vector2)
	{
		if (vector1.size() != vector2.size())
			return false;

		bool equal = false;
		unsigned int i = 0;

		do
		{
			equal = (vector1[i] == vector2[i]);
			++i;
		} while ((i < vector1.size()) && equal);
		
		return equal;
	}


	template <class T, int N> inline void output1OfN(const std::string &s)
	{
		static int i(0);
		if (!(i % N))
			std::cout << s << std::endl;		
		++i;
	}

	template <class T> void setupAll(std::vector<std::shared_ptr<T> > &v)
	{
		for (unsigned int i(0); i < v.size(); ++i)
			v[i]->setup();
	}

	template <class T> void computeAll(std::vector<std::shared_ptr<T> > &v)
	{
		for (unsigned int i(0); i < v.size(); ++i)
			v[i]->compute();
	}

	template <typename T1> class comparatorIndeces {
		private:
			std::vector<T1> & v1;
		public:	
			comparatorIndeces(std::vector<T1> & v): v1(v){}	
			inline bool operator()(unsigned int i, unsigned int j) {return v1[i]<v1[j];}
	};
	
	template <typename T1, typename T2> inline void sortTwoVectors(std::vector<T1> & v1, std::vector<T2> & v2)
	{
		if(v1.size()!=v2.size())
			errorMessage("sortTwoVectors: the vectors have different sizes");
		
		// creates index vector
		std::vector<unsigned int> ind(v1.size());
		for(unsigned int i(0); i<v1.size(); ++i) 
			ind[i]=i;
		
		sort(ind.begin(), ind.end(), comparatorIndeces<T1>(v1));
		reorderVector(ind,v1);
		reorderVector(ind,v2);
	}

	template <class T> inline void reorderVector(std::vector<unsigned int> & ind, std::vector<T> & v)
	{
		if(ind.size()!=v.size())
			errorMessage("reorderVector: the vectors have different sizes");

		unsigned int n(v.size());
		std::vector<T> nContainer(n);
		for(unsigned int i(0); i<n; ++i)
			nContainer[i]=v[ind[i]];
		v.swap(nContainer);
	}
	
	
} // asl namespace

#endif // ASLUTILITIES_H
