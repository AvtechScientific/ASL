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


#ifndef ACLMEMBLOCK_H
#define ACLMEMBLOCK_H

#include "../aclStdIncludes.h"
#include "../aclElementBase.h"
#include "../aclUtilities.h"
#include <aslUtilities.h>
#include "../aclHardware.h"


// using namespace asl;
// using shared_ptr;
namespace acl
{

	class MemBlock: public ElementBase
	{
		protected:
			shared_ptr<cl::Buffer> buffer; //< pointer in gpu
			weak_ptr<void> region;
			MemBlock();
			MemBlock(unsigned int size, TypeID typeID, CommandQueue queue_);
			MemBlock(unsigned int size, TypeID typeID, char *initArray, CommandQueue queue_);
			virtual void swapBuffers(MemBlock & a);
		public:
			virtual cl::Buffer &getBuffer();
			shared_ptr<void> map();
			friend inline void swapBuffers(MemBlock & a, MemBlock & b);
	};

	typedef shared_ptr<MemBlock> ElementData;
	inline void swapBuffers(MemBlock & a, MemBlock & b);

	
	template<typename T> inline std::shared_ptr<T> map(ElementData m);
	
// --------------------Implementation -----------------------

	inline void swapBuffers(MemBlock & a, MemBlock & b)
	{
		a.swapBuffers(b);
	}


	template<typename T> inline std::shared_ptr<T> map(ElementData m)
	{
		if (m->getTypeID() != typeToTypeID<T>())
			asl::errorMessage ("map: there is attempt to cast pointer with type " + 
			              typeToStr<T>() +
			              " for an element  with type " +
			              TYPE[m->getTypeID()]);
		return std::shared_ptr<T> (m->map(),(T*)m->map().get());
	}
} // namespace acl

#endif // ACLVECTOR_H
