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


#ifndef ASLBLOCKS_H
#define ASLBLOCKS_H

#include <iostream>
#include <fstream>
#include "../math/aslVectors.h"

namespace acl
{
	class VectorOfElements;
}

namespace asl
{

	/// The block without discretization (size and position are float)
	class Box
	{
		public:
			typedef AVec<> V; ///< Type of a vector
			V size;
			V position;

			/// the size and position are taken 0
			inline explicit Box(unsigned int nd);
	};
	
	
	/// Simple block with position
	/// \ingroup DataFields 
	/// The Block describes the simulation grid
	/// and is defined by the position and by the size.
	class Block
	{
		public:
			typedef AVec<int> DV; ///<Discrete Vector
		private:	
			DV size;		
		public:
			typedef AVec<> V; ///<Type of the position

			V position;
			double dx;
			DV c2iTransformVector;

			/// the size is taken 1, the position is taken to be 0
			inline Block();
			inline explicit Block(unsigned int nd);
			inline Block(const DV & s, double dx, const V & p);
			inline explicit Block(const DV & s, double dx = 1);
			inline Block(const Block & b);	
			inline const Block& operator=(const Block & b); 
			/// defines convertion rule of 1D/2D/3D index \p i into container one
			int c2i(const Block::DV & c) const;
			inline void setSize(const DV & s);
			inline const DV & getSize() const;
			/// returns position of the point diagonal to (0) one
			/// returns Body Diagonal (or space diagonal) that starts at (0)
			inline const V getBPosition() const;
			acl::VectorOfElements initACLPositionDiscrete();
			acl::VectorOfElements initACLPosition();

			acl::VectorOfElements getACLPositionDiscrete();
			acl::VectorOfElements getACLPosition();
	};


	/// results Block which is inset or outset of the \p bl on value of a
	/// \ingroup DataFields
	const Block offset(const Block & bl, int a = 1);


	/// Dimensionality of the block \ingroup DataFields
	/// \todo rename here and everywhere to dimNum???
	inline const unsigned int nD(const Block & b);

	/// Checks whether \p a inside of \p b
	/// \ingroup DataFields
	inline const bool in(const Block & b, AVec<> a);

	/// Checks whether \p a inside of \p b
	/// \ingroup DataFields
	inline const bool in(const Block & b, AVec<int> a);

	
	inline const AVec<int> continiousToDiscret(const Block & b, AVec<> a);

	//--------------- Implementation----------------

	inline Box::Box(unsigned int nd):size(nd),position(nd)
	{
	}	


	inline AVec<int> castTransformVector(AVec<int> s)
	{
		unsigned int n(s.getSize());
		AVec<int> r(n, 1);
		int a(1);
		for (unsigned int i = 0; i < n - 1; ++i)
		{
			a *= s[n - 1 - i];
			r[n - 2 - i] = a;			
		}
		return r;
	}
	

	inline Block::Block(unsigned int nd):
		size(nd),
		position(nd),
		dx(1),
		c2iTransformVector(castTransformVector(size))
	{
	}


	inline Block::Block():
		size(1),
		position(1,0.),
		dx(1),
		c2iTransformVector(castTransformVector(size))
	{
	}
	

	inline Block::Block(const Block::DV & s, double d):
		size(s),
		position(V(s.getSize())),
		dx(d),
		c2iTransformVector(castTransformVector(s))
	{
	}


	inline Block::Block(const Block::DV & s, double d, const Block::V & p):
		size(s),
		position(p),
		dx(d),
		c2iTransformVector(castTransformVector(s))			
	{
		if (p.getSize() != s.getSize())
			errorMessage ("Block::Block() Size and Position dimensionalities are different");
	}


	inline Block::Block(const Block & b):
		size(b.size),
		position(b.position),
		dx(b.dx),
		c2iTransformVector(b.c2iTransformVector)		
	{
	}


	inline const Block & Block::operator=(const Block & b)
	{
		position=b.position;
		size = b.size;
		dx = b.dx;
		c2iTransformVector = b.c2iTransformVector;
		return *this;
	}


	inline int Block::c2i(const DV & c) const
	{
		unsigned int n(c2iTransformVector.getSize());
		if (c.getSize() != n) 
			errorMessage("Block::c2i() - The input vector size does not correspond to the block dimensionality");
		
		return c * c2iTransformVector;
	}
		

	inline void Block::setSize(const DV & s)
	{
		size = s;
		if(position.getSize()!= s.getSize())
			position=V(s.getSize());
		c2iTransformVector = castTransformVector(s);
	}


	inline const Block::DV & Block::getSize() const
	{
		return size;
	}


	inline const Block::V Block::getBPosition() const
	{
		return position + (V(size - DV(size.getSize(), 1.)) ) * dx;
	}


	inline const unsigned int nD(const Block & b)
	{
		return b.position.getSize();
	}
		

	inline const bool in(const Block & b, AVec<> a)
	{
		return positive(a-b.position) && positive(b.getBPosition() - a);   
	}


	inline const bool in(const Block & b, AVec<int> a)
	{
		return positive(a) && positive(b.getSize() - a);   
	}
				
}	//asl

#endif //ASLBLOCKS_H

