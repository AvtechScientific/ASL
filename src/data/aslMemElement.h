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


#ifndef ASLMEMELEMENT
#define ASLMEMELEMENT

#include <memory>


using  namespace std;

namespace asl
{

	class MemElementBase
	{
		protected:
			unsigned int size;
			
			inline MemElementBase(unsigned int n = 0);
		public:			
			virtual ~MemElementBase() = 0;
			inline unsigned int getSize() const;
			virtual void resize() = 0;
			
	};

	typedef std::shared_ptr<MemElementBase> MemElement;

	template <typename T> class MemVector: public MemElementBase
	{
		private:
			T* container;
			bool createdContainer;
		public:
			MemVector(unsigned int n = 0);
			virtual ~MemVector();			
			virtual void resize(unsigned int n);
			void setContainer(unsigned int n, T* p);
	};


//------------------------ Implementations -----------------------

	inline MemElementBase::MemElementBase(unsigned int n): size(n)
	{
	}

	inline unsigned int MemElementBase::getSize() const
	{
		return size;
	}
	


} //asl
#endif //ASLMEMELEMNT

