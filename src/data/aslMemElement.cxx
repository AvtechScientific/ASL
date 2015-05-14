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


#include "aslMemElement.h"

namespace asl{

	template <typename T> MemVector<T>::MemVector(unsigned int n):
		MemElementBase(n),container(new T[n]),createdContainer(true)
	{
	}

	template <typename T> MemVector<T>::~MemVector()
	{
		if (createdContainer)
			delete[] container;
	}

	template <typename T> void MemVector<T>::resize(unsigned int n)
	{
		if (createdContainer)
			delete[] container;
		container=new T[n];
		createdContainer=true;
		size=n;
	}

	template <typename T> void MemVector<T>::setContainer(unsigned int n, T* p)
	{
		if (createdContainer)
			delete[] container;
		container=p;
		createdContainer=false;
		size=n;
	}
	

}//asl


