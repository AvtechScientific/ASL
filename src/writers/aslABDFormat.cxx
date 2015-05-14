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


#include "aslABDFormat.h"

#include "data/aslBlocks.h"
#include <acl/acl.h>
#include <utilities/aslSmartPtrUtils.h>
#include <acl/DataTypes/aclArray.h>

using  namespace std;

namespace asl {

	ABDFileOut & operator <<(ABDFileOut & f, const Block &b)
	{
		f<<b.dx<<b.getSize()<<b.position;		
		return f;
	}

	ABDFileIn & operator >>(ABDFileIn & f,Block &b)
	{
		double dx(0);
		Block::V pos;
		Block::DV size;
		f>>dx>>size>>pos;

		b.dx=dx;
		b.position=pos;
		b.setSize(size);
		
		return f;
	}

	ABDFileOut & operator <<(ABDFileOut & f, const AbstractData &a)
	{
		unsigned int s(a.getDContainer()[0]->getSize());
		double* d=new double[s];
		copy(a.getDContainer()[0],d);
		f<<std::make_pair(d,s);
		delete[] d;
		return f;
	}

	ABDFileIn & operator >>(ABDFileIn & f, AbstractData &a)
	{
		std::shared_ptr<double> d(acl::map<double>(a.getDContainer()[0]));
		unsigned int s(a.getDContainer()[0]->getSize());
		f>>std::make_pair(d.get(),s);
		return f;
	}

	ABDFileIn & get(ABDFileIn & f, AbstractData &a, std::shared_ptr<double> d)
	{
		unsigned int s(a.getDContainer()[0]->getSize());
		if (d.get()==0)
			d.reset(new double[s],ArrayDeleter<double>());
		f>>std::make_pair(d.get(),s);
		copy(d.get(),a.getDContainer()[0]);
		return f;
	}
	
	
} //asl


