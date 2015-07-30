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


#ifndef ACLSUBVECTOR_H
#define ACLSUBVECTOR_H

//#include <CL/cl.hpp>
// Supply "cl.hpp" with ASL, since it is not present in OpenCL 2.0
// Remove the file after switching to OpenCL 2.1
#include "acl/cl.hpp"
#include "aclArray.h"
#include "aslUtilities.h"

///\todo{Add padding at the end?}

namespace acl
{

	template <typename T> class Subvector: public MemBlock
	{
		private:
			string name;
			static const string prefix;
			static unsigned int id;			
			unsigned int offset;
			shared_ptr<Array<T> > vector;
			cl_buffer_region buffer_create_info;
		public:
			Subvector(shared_ptr<Array<T> > vector_,
			          unsigned int size_,
			          unsigned int offset_);
			virtual cl::Buffer &getBuffer();
			virtual string str(const KernelConfiguration & kernelConfig) const;
			virtual string getName() const;
			virtual string getTypeSignature(const KernelConfiguration & kernelConfig) const;
			virtual string getLocalDeclaration(const KernelConfiguration & kernelConfig) const;
			virtual void addToKernelSource(std::vector<Element> & arguments,
			                               std::vector<Element> & localDeclarations) const;
			virtual void setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const;			
	};



//---------------------------- Implementations------------------------

	template <typename T> Subvector<T>::Subvector(shared_ptr<Array<T> > vector_,
	                                              unsigned int size_,
	                                              unsigned int offset_):
		MemBlock(),
		offset(offset_),
		vector(vector_)
	{
		size = size_;
		queue = vector_->getQueue();
		if ( (offset + size) > vector->getSize() )
		{
			errorMessage("Subvector::Subvector() - (offset + size) > vector->getSize()");
		}
		else
		{
			buffer_create_info.origin = offset * sizeof(T);
			buffer_create_info.size = size * sizeof(T);
		}
		
		++id;
		name = prefix + asl::numToStr(id);
	}


	template <typename T> cl::Buffer &Subvector<T>::getBuffer()
	{
		cl_int status = 0;
		///  if original buffer is moved within GPU memory	
		*buffer = vector->getBuffer().createSubBuffer(CL_MEM_READ_WRITE,
		                                              CL_BUFFER_CREATE_TYPE_REGION,
		                                              &buffer_create_info,
		                                              &status);
		errorMessage(status, "Subvector::Subvector() - createSubBuffer()");

		return *buffer;
	}


	template <typename T> string Subvector<T>::getName() const
	{
		return name;
	}


	template <typename T> string Subvector<T>::str(const KernelConfiguration & kernelConfig) const
	{
		return name + "[" + INDEX + "]";
	}


	template <typename T> string Subvector<T>::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "__global " + typeToStr<T>(kernelConfig.vectorWidth) + " *" + name;
	}


	template <typename T> string Subvector<T>::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	// Must be empty. Only operators can add arguments.
	template <typename T> 
	void Subvector<T>::addToKernelSource(std::vector<Element> & arguments,
	                                     std::vector<Element> & localDeclarations) const
	{
	}


	template <typename T> 
	void Subvector<T>::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
	{
		cl_int status = 0;		
		status = kernel.setArg(argumentIndex, *buffer);
		errorMessage(status, "Kernel::setArg() - " + name + ", argument " + numToStr(argumentIndex));
	}

} // namespace acl

#endif // ACLSUBVECTOR_H
