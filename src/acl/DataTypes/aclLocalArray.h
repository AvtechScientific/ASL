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


#ifndef ACLLOCALARRAY_H
#define ACLLOCALARRAY_H

#include "../aclElementBase.h"
#include "../../aslUtilities.h"
#include "../aclUtilities.h"

using namespace asl;

namespace acl
{


	/** 
		This Element appears in the Kernel source code as a local array.
		It is not passed to the Kernel as an argument.
	*/
	template <typename T> class LocalArray: public ElementBase
	{
		private:
			string name;
			static const string prefix;
			static unsigned int id;
		public:
			explicit LocalArray(unsigned int size_);
			virtual string str(const KernelConfiguration & kernelConfig) const;
			virtual string getName() const;
			virtual string getAddressSpaceQualifier() const;
			virtual string getTypeSignature(const KernelConfiguration & kernelConfig) const;
			virtual string getLocalDeclaration(const KernelConfiguration & kernelConfig) const;
			virtual void addToKernelSource(vector<Element> & arguments,
			                               vector<Element> & localDeclarations) const;
			virtual void setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const;
	};


	template <typename T> LocalArray<T>::LocalArray(unsigned int size_):
		ElementBase(true, size_, typeToTypeID<T>())
	{
		++id;
		name = prefix + asl::numToStr(id);
	}


	template <typename T> string LocalArray<T>::str(const KernelConfiguration & kernelConfig) const
	{
		if (kernelConfig.unaligned && kernelConfig.vectorWidth > 1)
		{
			return "vload" + numToStr(kernelConfig.vectorWidth) + "(0, &" + name + "[" + INDEX + "])";
		}
		else
		{
			return name + "[" + INDEX + "]";
		}
	}


	template <typename T> string LocalArray<T>::getName() const
	{
		return name;
	}
		

	template <typename T> string LocalArray<T>::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	template <typename T> string LocalArray<T>::getAddressSpaceQualifier() const
	{
		return "__local";
	}


	template <typename T> string LocalArray<T>::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		return "__local " + typeToStr<T>(kernelConfig.unaligned ? 1 : kernelConfig.vectorWidth)
				+ " " + name + "[" + asl::numToStr(size + paddingElements(size, kernelConfig)) + "]";
	}


	template <typename T> void LocalArray<T>::addToKernelSource(vector<Element> & arguments,
	                                                             vector<Element> & localDeclarations) const
	{
	}


	// can't be set as an argument
	template <typename T> void LocalArray<T>::setAsArgument(cl::Kernel & kernel,
	                                                         unsigned int argumentIndex) const
	{
	}


} // namespace acl

#endif // ACLLOCALARRAY_H
