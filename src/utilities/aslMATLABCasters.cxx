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


#include "aslMATLABCasters.h"

#include <data/aslBlocks.h>
#include <data/aslProbe.h>
#include <data/aslDataWrapper.h>
#include <acl/acl.h>
#include <acl/DataTypes/aclArray.h>
#include <math/aslVectors.h>

#include <matio.h>

template <typename T> asl::AVec<T> castMATLABVector(const asl::AVec<T> &a)
{
	unsigned int n(nD(a));
	asl::AVec<T> b(n>2?n:2);
	if(n == 1)
		b[0] = a[0]; b[1] = 1;
	if (n >= 2)
		b = a;
	return b;
}

template <typename T> class MATLABTypes{};
template <> class MATLABTypes<double>{public: 
		const static int c=MAT_C_DOUBLE,
						 t=MAT_T_DOUBLE;
};
template <> class MATLABTypes<float>{public: 
		const static int c=MAT_C_SINGLE,
						 t=MAT_T_SINGLE;
};
template <> class MATLABTypes<cl_int>{public: 
		const static int c=MAT_C_INT32,
						 t=MAT_C_INT32;
};



matiocpp::Var::~Var(){
	if (var->mem_conserve && freeArray)
		switch (var->data_type)
		{
			case MAT_T_DOUBLE: delete[] ((double*)var->data); break;
			case MAT_T_SINGLE: delete[] ((float*)var->data); break;
			case MAT_T_INT32: delete[] ((cl_int*)var->data); break;
		}
	Mat_VarFree(var);
}

template<typename T> matvar_t* matiocpp::Var::init(T *d, const asl::AVec<int> & size, const std::string &name)
{
	asl::AVec<size_t> msize(castMATLABVector(size));
	return Mat_VarCreate(name.c_str(),
	                 	 MATLABTypes<T>::c,
	                   	 MATLABTypes<T>::t,
	                   	 nD(msize),
	                   	 &msize[0],
	                   	 d,
	                   	 0);
}

template<typename T> matiocpp::Var::Var(T *d, const asl::AVec<int> & size, const std::string &name, bool freeArr):
	var(init(d,size,name)), 
	freeArray(freeArr)
{
}

template<typename T> matiocpp::Var::Var(T *d, unsigned int size, const std::string &name, bool freeArr):
	var(init(d,asl::makeAVec((int)size),name)),
	freeArray(freeArr)
{	
}

template<typename T> matiocpp::Var::Var(std::vector<T> & d, const std::string &name):
	var(init(&d.front(),asl::makeAVec((int)d.size()),name)),
	freeArray(false)
{	
}


namespace asl
{

	SPMatVar castMATLABCellArray(acl::Element source, const std::string &name)
	{
		if (!isMemBlock(source))
			errorMessage ("Error (castMATLABCellArray): the given element is not a MemBlock type");
		
		if (source->getTypeID() == acl::TYPE_DOUBLE)
		{
			double* d=new double[source->getSize()];
			copy(source,d);
			return SPMatVar(new matiocpp::Var(d,source->getSize(),name));
		}
		if (source->getTypeID() == acl::TYPE_FLOAT)
		{
			float* d=new float[source->getSize()];
			copy(source,d);
			return SPMatVar(new matiocpp::Var(d,source->getSize(),name));
		}
		if (source->getTypeID() == acl::TYPE_INT)
		{
			int* d=new int[source->getSize()];
			copy(source,d);
			return SPMatVar(new matiocpp::Var(d,source->getSize(),name));
		}
		return SPMatVar();
	}	

	SPMatVar castMATLABCellArray(acl::Element source, const AVec<int> & size, const std::string &name)
	{
		if (!isMemBlock(source))
			errorMessage ("Error (castMATLABCellArray): the given element is not a MemBlock type");
		
		if (source->getTypeID() == acl::TYPE_DOUBLE)
		{
			double* d = new double[source->getSize()];
			copy(source, d);
			return SPMatVar(new matiocpp::Var(d,size,name));
		}
		if (source->getTypeID() == acl::TYPE_FLOAT)
		{
			float* d = new float[source->getSize()];
			copy(source, d);
			return SPMatVar(new matiocpp::Var(d,size,name));
		}
		if (source->getTypeID() == acl::TYPE_INT)
		{
			int* d = new int[source->getSize()];
			copy(source, d);
			return SPMatVar(new matiocpp::Var(d,size,name));
		}
		return SPMatVar();
	}	


	SPMatVar castMATLABCellArray(const AbstractData & d, const std::vector<std::string> &names)
	{
		const Block &b(d.getBlock());
		bool bnames(names.size()>0);

		if (d.getDContainer().size()==1)
			return castMATLABCellArray(d.getDContainer()[0],b.getSize(),bnames?names[0]:"");
//		else
//			for (unsigned int i(1); i < ; ++i)
//			image->GetPointData()->AddArray(castVTKDataArray(d.getDContainer()[i],bnames?names[i]:""));
		return SPMatVar();				
	}

	SPMatVar castMATLABCellArray(Probe &p, unsigned int component, const std::string &name)
	{
		return SPMatVar(new matiocpp::Var(p.getComponent(component),name));
	}	
	

}// asl
