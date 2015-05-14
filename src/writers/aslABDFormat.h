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


#ifndef ASLABDFORMAT_H
#define ASLABDFORMAT_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <string>
#include <memory>

#include <data/aslDataWrapper.h>

using  namespace std;

namespace asl {

	/**
		\defgroup ABDFormat ASL Binary Dump (ABD) format
		\ingroup IO
	 */
	
	/// ABD (ASL Binary Dump) file, input
	/**
		 \ingroup ABDFormat
	*/
	class ABDFileIn: public std::ifstream{
		public:
			inline ABDFileIn();
			inline ABDFileIn(string name);
	};

	/// ABD (ASL Binary Dump) file, output
	/**
		 \ingroup ABDFormat
	*/
	class ABDFileOut: public std::ofstream{
		public:
			inline ABDFileOut();
			inline ABDFileOut(string name);
	};
	
/*	/// \relates ABDFileOut
	template <typename T> inline ABDFileOut & operator <<(ABDFileOut & f, const T & a);
	/// \relates ABDFileIn
	template <typename T> inline ABDFileIn & operator >>(ABDFileIn & f,T & a);
*/
	
	/// \relates ABDFileOut
	inline ABDFileOut & operator <<(ABDFileOut & f, const int a);
	/// \relates ABDFileIn
	inline ABDFileIn & operator >>(ABDFileIn & f,int & a);
	/// \relates ABDFileOut
	inline ABDFileOut & operator <<(ABDFileOut & f, const unsigned int a);
	/// \relates ABDFileIn
	inline ABDFileIn & operator >>(ABDFileIn & f, unsigned int & a);
	/// \relates ABDFileOut
	inline ABDFileOut & operator <<(ABDFileOut & f, const float a);
	/// \relates ABDFileIn
	inline ABDFileIn & operator >>(ABDFileIn & f,float & a);
	/// \relates ABDFileOut
	inline ABDFileOut & operator <<(ABDFileOut & f, const double a);
	/// \relates ABDFileIn
	inline ABDFileIn & operator >>(ABDFileIn & f,double & a);
	
	
	/// \relates ABDFileOut
	template <typename T> inline ABDFileOut & operator <<(ABDFileOut & f,pair<T*, unsigned int> a);
	/// \relates ABDFileIn
	template <typename T> inline ABDFileIn & operator >>(ABDFileIn & f,pair<T*, unsigned int> a);

	/// \relates ABDFileOut	
	inline ABDFileOut & operator <<(ABDFileOut & f, const string &a);
	/// \relates ABDFileIn	
	inline ABDFileIn & operator >>(ABDFileIn & f,string &a);

	/// \relates ABDFileOut
	template <typename T> inline ABDFileOut & operator <<(ABDFileOut & f, const AVec<T> & a);
	/// \relates ABDFileIn
	template <typename T> inline ABDFileIn & operator >>(ABDFileIn & f,AVec<T> & a);

	class Block;
	
	/// \relates ABDFileOut	
	ABDFileOut & operator <<(ABDFileOut & f, const Block &a);
	/// \relates ABDFileIn	
	ABDFileIn & operator >>(ABDFileIn & f, Block &a);

	class AbstractData;
	
	/// writes data. It is assumed that the Block is written separately \relates ABDFileOut	
	/**
		 The function writes only the first element
	*/
	ABDFileOut & operator <<(ABDFileOut & f, const AbstractData &a);

	/// reads data. It is assumed that the class has the propper size \relates ABDFileIn	
	/**
		 The function writes only the first element
	*/
	ABDFileIn & operator >>(ABDFileIn & f, AbstractData &a);

	/// reads data. It is assumed that the class has the propper size \relates ABDFileIn	
	/**
		 The function reads only the first element. Additionaly it creates and 
		 stores data in the memory \p d; in case of d.get()==0 defineds whether 
		 the new d should be allocated
	*/
	ABDFileIn & get(ABDFileIn & f, AbstractData &a, std::shared_ptr<double> d);
		 		
	/// writes \p data in a file with ABD (ASL Binary Dump) format
	/**
		 \ingroup ABDFormat
	*/
	void writeABD(const string &fileName, const AbstractData & data, const string & name);


//------------------------IMPLEMENTATION------------------------------

	ABDFileIn::ABDFileIn():ifstream()
	{};
	ABDFileIn::ABDFileIn(string name):
		ifstream(name,ios::in | ios::binary)
	{}	

	ABDFileOut::ABDFileOut():ofstream()
	{}
	ABDFileOut::ABDFileOut(string name):
		ofstream(name,ios::out | ios::binary)
	{}

	inline ABDFileOut & operator <<(ABDFileOut & f, const int a)
	{
		f.write((char*)&a,sizeof(int)); 
		return f;
	}

	inline ABDFileIn & operator >>(ABDFileIn & f,int & a)
	{
		f.read((char*)&a,sizeof(int)); 
		return f;
	}

	inline ABDFileOut & operator <<(ABDFileOut & f, const unsigned int a)
	{
		f.write((char*)&a,sizeof(unsigned int)); 
		return f;
	}

	inline ABDFileIn & operator >>(ABDFileIn & f, unsigned int & a)
	{
		f.read((char*)&a,sizeof(unsigned int)); 
		return f;
	}

	inline ABDFileOut & operator <<(ABDFileOut & f, const float a)
	{
		f.write((char*)&a,sizeof(float)); 
		return f;
	}

	inline ABDFileIn & operator >>(ABDFileIn & f,float & a)
	{
		f.read((char*)&a,sizeof(float)); 
		return f;
	}

	inline ABDFileOut & operator <<(ABDFileOut & f, const double a)
	{
		f.write((char*)&a,sizeof(double)); 
		return f;
	}

	inline ABDFileIn & operator >>(ABDFileIn & f,double & a)
	{
		f.read((char*)&a,sizeof(double)); 
		return f;
	}
		

	template <typename T> inline ABDFileOut & operator <<(ABDFileOut & f,pair<T*, unsigned int> a)
	{
		f.write((char*)a.first,sizeof(T)*a.second); 
		return f;
	}
	template <typename T> inline ABDFileIn & operator >>(ABDFileIn & f,pair<T*, unsigned int> a)
	{
		f.read((char*)a.first,sizeof(T)*a.second); 
		return f;
	}
	
	inline ABDFileOut & operator <<(ABDFileOut & f, const string &a){
		unsigned int n=a.size();
		f<<n<<make_pair(&a[0],n);
		return f;
	}

	inline ABDFileIn & operator >>(ABDFileIn & f,string &a){
		unsigned int n;
		f>>n; a.resize(n);
		f>>make_pair(a.data(),n);
		return f;
	}

	template <typename T> inline ABDFileOut & operator <<(ABDFileOut & f, const AVec<T> & a)
	{
		unsigned int n(a.getSize());
		f<<n<<make_pair(&(a[0]),n);
		return f;
	}

	template <typename T> inline ABDFileIn & operator >>(ABDFileIn & f,AVec<T> & a)
	{
		unsigned int n(0);
		f>>n; 
		a.resize(n);
		f>>make_pair(&(a[0]),n);
		return f;
	}

		
}// asl

#endif //ASLVTKFORMAT_H

