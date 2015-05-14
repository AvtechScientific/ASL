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


#ifndef ACLARRAY_H
#define ACLARRAY_H

#include "../aclStdIncludes.h"
#include "aclMemBlock.h"

using namespace asl;

namespace acl
{

	/// Global array
	template <typename T> class Array: public MemBlock
	{
		protected:
			string name;
			static const string prefix;
			static unsigned int id;
		public:
			explicit Array(unsigned int size, CommandQueue queue_ = hardware.defaultQueue);
			Array(unsigned int size, T *initArray, CommandQueue queue_ = hardware.defaultQueue);
			virtual string str(const KernelConfiguration & kernelConfig) const;
			virtual string getName() const;
			virtual string getAddressSpaceQualifier() const;
			virtual string getTypeSignature(const KernelConfiguration & kernelConfig) const;
			virtual string getLocalDeclaration(const KernelConfiguration & kernelConfig) const;
			virtual void addToKernelSource(vector<Element> & arguments,
			                               vector<Element> & localDeclarations) const;
			virtual void setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const;
			/// swaps buffer values of this vector and vector \p a;
	};

	typedef shared_ptr<Array<cl_int> > ElementArrayInt;
	typedef shared_ptr<Array<cl_float> > ElementArrayFloat;
	typedef shared_ptr<Array<cl_double> > ElementArrayDouble;
	typedef shared_ptr<Array<cl_long> > ElementArrayLong;
	
//------------------- Implementation -------------------------------
	
	template <typename T> Array<T>::Array(unsigned int size_, CommandQueue queue_):
		MemBlock(size_, typeToTypeID<T>(), queue_)
	{
		++id;
		name = prefix + asl::numToStr(id);
	}


	template <typename T> Array<T>::Array(unsigned int size_, T *initArray, CommandQueue queue_):
		MemBlock(size_, typeToTypeID<T>(), queue_, (char*)initArray)
	{
		++id;
		name = prefix + asl::numToStr(id);
	}



	template <typename T> string Array<T>::getName() const
	{
		return name;
	}


	template <typename T> string Array<T>::str(const KernelConfiguration & kernelConfig) const
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


	template <typename T> string Array<T>::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "__global " + typeToStr<T>(kernelConfig.unaligned ? 1 : kernelConfig.vectorWidth) + " *" + name;
	}


	template <typename T> string Array<T>::getAddressSpaceQualifier() const
	{
		return "__global";
	}
	

	template <typename T> string Array<T>::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	// Must be empty. Only operators can add arguments.
	template <typename T> 
	void Array<T>::addToKernelSource(vector<Element> & arguments, vector<Element> & localDeclarations) const
	{
	}


	template <typename T> 
	void Array<T>::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
	{
		cl_int status = 0;		
		status = kernel.setArg(argumentIndex, *buffer);
		errorMessage(status, "Kernel::setArg() - " + name
		             + ", argument " + numToStr(argumentIndex));
	}


} // namespace acl

#endif // ACLARRAY_H
