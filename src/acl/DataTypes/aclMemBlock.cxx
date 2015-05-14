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


#include "aclMemBlock.h"
#include "../../aslUtilities.h"
#include "../aclUtilities.h"


using namespace asl;

namespace acl
{

	/// class for compatibility with std::shared_ptr enssure correct unmapping \ingroup LDI
	/// \todoIgal: maybe pass a MemBlock to the constructor instead of a buffer and a queue?
	template <typename T> class VectorUnmapper
	{
			shared_ptr<cl::Buffer> buffer;
			CommandQueue queue;
		public:
			VectorUnmapper(shared_ptr<cl::Buffer> buffer_, CommandQueue queue_) :
				buffer(buffer_),
				queue(queue_)
			{};
			
			inline void operator() (T* d) const
			{
				cl_int status = 0;
				cl::Event event;
				
				status = queue->enqueueUnmapMemObject(*buffer,
						                              d,
						                              NULL,
						                              &event);
				errorMessage(status, "enqueueUnmapMemObject()");
				status = event.wait();
				errorMessage(status, "Event::wait() - event");			
				d = NULL;
			}
	};


	MemBlock::MemBlock():
		// typeID fictitious
		ElementBase(true, 0, TYPE_INT)
	{
	    buffer.reset(new cl::Buffer());
	}


	MemBlock::MemBlock(unsigned int size_, TypeID typeID, CommandQueue queue_):
		ElementBase(true, size_, typeID)
	{
		queue = queue_;

		cl_int status = 0;

		// allocate value in GPU memory with padding
	    buffer.reset(new cl::Buffer(getContext(queue),
	                                CL_MEM_READ_WRITE,
	                                size * TYPE_SIZE[typeID] + paddingBytes(size, TYPE_SIZE[typeID], queue),
	                                NULL,
	                                &status));
		errorMessage(status, "cl::Buffer()");
	}


	MemBlock::MemBlock(unsigned int size_, TypeID typeID, char *initArray, CommandQueue queue_):
		ElementBase(true, size_, typeID)
	{
		queue = queue_;
		
		cl_int status = 0;

		// allocate value in GPU memory with padding
	    buffer = shared_ptr<cl::Buffer>(new cl::Buffer(getContext(queue),
	                                                   CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
	                                                   size * TYPE_SIZE[typeID] + paddingBytes(size, TYPE_SIZE[typeID], queue),
	                                                   initArray,
	                                                   &status));
		errorMessage(status, "cl::Buffer()");
	}


	cl::Buffer & MemBlock::getBuffer()
	{
		return *buffer;
	}


	std::shared_ptr<void> MemBlock::map()
	{
		shared_ptr<void> p;
		if (region.expired())
		{
			cl_int status = 0;
			cl::Event event;

			// blocking_map = CL_TRUE (second argument)
			// read and write access (third argument)
			p.reset(queue->enqueueMapBuffer(*buffer,
				                            CL_TRUE,
				                            CL_MAP_READ | CL_MAP_WRITE,
				                            0,
				                            size * TYPE_SIZE[typeID],
				                            NULL,
				                            &event,
				                            &status),
			         VectorUnmapper<void> (buffer, queue));
			errorMessage(status, "enqueueMapBuffer()");
			status = event.wait();
			errorMessage(status, "Event::wait() - event");			
			region = p;
		}
		return region.lock();
	}


	void MemBlock::swapBuffers(MemBlock & a)
	{
		if (compatible(size, queue, a.getSize(), a.getQueue()))
		{
			std::swap(buffer, a.buffer);
		}
		else
		{
			errorMessage("acl::Vector::swap() - operands are incompatible. \
				 		 Either they reside on different devices or their sizes do not match: "
			             + numToStr(size) + " and " + numToStr(a.size));
		}
	}
}
