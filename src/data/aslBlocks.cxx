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


#include "aslBlocks.h"

namespace asl{
	
	const Block offset(const Block & bl, int a)
	{
		if (minComponent(bl.getSize())<-a*2) 
			errorMessage("The inset is larger than the block size");
		Block b(bl.getSize().getSize());
		b.position=bl.position-Block::V(bl.getSize().getSize(),a*bl.dx);
		b.setSize(bl.getSize()+Block::DV(b.getSize().getSize(),a*2));
		b.dx=bl.dx;
		return b;
	}

}// asl
