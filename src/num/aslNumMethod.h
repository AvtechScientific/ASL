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


#ifndef ASLNUMMETHOD_H
#define ASLNUMMETHOD_H

#include <acl/aclStdIncludes.h>
#include <memory>

namespace asl
{
	/// Virtual class describes general interface for Numerical methods;
	/// \ingroup Numerics
	class NumMethod
	{
		public:
			/// Executes the numerical procedure
			virtual void execute() = 0;
			/// Builds the necesery internal data and kernels
			virtual void init() = 0;
			virtual ~NumMethod();
	};

	typedef std::shared_ptr<NumMethod> SPNumMethod;


	template <class T> inline void initAll(std::vector<T*> &v);
	template <class T> inline void initAll(std::vector<std::shared_ptr<T> > & v);
	
	template <class T> inline void executeAll(std::vector<T*> &v);
	template <class T> inline void executeAll(std::vector<std::shared_ptr<T> > & v);

//---------------------------- Implementation ---------------------------

	template <class T> void executeAll(std::vector<T*> &v)
	{
		for (unsigned int i(0); i < v.size(); ++i)
			v[i]->execute();
	}

	template <class T> void executeAll(std::vector<std::shared_ptr<T> > &v)
	{
		for (unsigned int i(0); i < v.size(); ++i)
			v[i]->execute();
	}

	template <class T> void initAll(std::vector<T*> &v)
	{
		for (unsigned int i(0); i < v.size(); ++i)
			v[i]->init();
	}

	template <class T> void initAll(std::vector<std::shared_ptr<T> > &v)
	{
		for (unsigned int i(0); i < v.size(); ++i)
			v[i]->init();
	}

	
}	//asl

#endif //ASLNUMMETHOD_H
