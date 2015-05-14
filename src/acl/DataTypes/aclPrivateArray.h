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


#ifndef ACLPRIVATEARRAY_H
#define ACLPRIVATEARRAY_H

#include "../aclElementBase.h"
#include "../../aslUtilities.h"
#include "../aclUtilities.h"

using namespace asl;

namespace acl
{


	/** 
		This Element appears in the Kernel source code
	 	as a private array (without address space qualifier).
		It is not passed to the Kernel as an argument.
	*/
	template <typename T> class PrivateArray: public ElementBase
	{
		private:
			string name;
			static const string prefix;
			static unsigned int id;
			vector<T> initVector;
		public:
			explicit PrivateArray(const vector<T> & initVector_);
			virtual string str(const KernelConfiguration & kernelConfig) const;
			virtual string getName() const;
			virtual string getAddressSpaceQualifier() const;
			virtual string getTypeSignature(const KernelConfiguration & kernelConfig) const;
			virtual string getLocalDeclaration(const KernelConfiguration & kernelConfig) const;
			virtual void addToKernelSource(vector<Element> & arguments,
			                               vector<Element> & localDeclarations) const;
			virtual void setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const;
	};


	template <typename T> PrivateArray<T>::PrivateArray(const vector<T> & initVector_):
		ElementBase(true, initVector_.size(), typeToTypeID<T>()),
		initVector(initVector_)
	{
		++id;
		name = prefix + asl::numToStr(id);
	}


	template <typename T> string PrivateArray<T>::str(const KernelConfiguration & kernelConfig) const
	{
		if (kernelConfig.vectorWidth > 1)
		{
			errorMessage("PrivateArray should not be used in a SIMD Kernel");
			return "";
		}
		else
			return name + "[" + INDEX + "]";
	}


	template <typename T> string PrivateArray<T>::getName() const
	{
		return name;
	}
		

	template <typename T> string PrivateArray<T>::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	template <typename T> string PrivateArray<T>::getAddressSpaceQualifier() const
	{
		return "__private";
	}


	template <typename T> string PrivateArray<T>::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		string s = typeToStr<T>() + " " + name + "[" + asl::numToStr(size) + "] = {";

		for (unsigned int i = 0; i < size; ++i)
			s += numToStr(initVector[i]) + ", ";

		// drop last ", "
		s.erase(s.size() - 2);
		s += "}";
			
		return s;
	}


	template <typename T> void PrivateArray<T>::addToKernelSource(vector<Element> & arguments,
	                                                              	                                                              vector<Element> & localDeclarations) const
	{
	}


	// can't be set as an argument
	template <typename T> void PrivateArray<T>::setAsArgument(cl::Kernel & kernel,
	                                                         unsigned int argumentIndex) const
	{
	}


} // namespace acl

#endif // ACLPRIVATEARRAY_H
