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


#include "aslVTKCasters.h"

#include <data/aslBlocks.h>
#include <data/aslDataWrapper.h>
#include <acl/acl.h>
#include <acl/DataTypes/aclArray.h>
#include <data/aslDataWithGhostNodes.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkLongLongArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>

#include <vtkImageData.h>

#include "acl/aclTypesList.h"


namespace asl
{

	template <typename T> void combineArrays(T* d1, unsigned int size, T* dTarget, unsigned int nComponents=1);
	template <typename T> void combineArrays(T* d1, T* d2, unsigned int size, T* dTarget,unsigned int nComponents=2);	
	template <typename T> void combineArrays(T* d1, T* d2, T* d3, unsigned int size, T* dTarget,unsigned int nComponents=3);
	template <typename TO, typename TI> TO* combineArrays(TI* d0, TI* d1, TI* d2, TI* d3, unsigned int size);

	template <typename T> void combineArrays(T* d1,
	                                         unsigned int size,
	                                         T* dTarget,
	                                         unsigned int nComponents = 1)
	{
		if (nComponents < 1)
			errorMessage("combineArrays() - attempt to provide nComponents that is less than 1");

		for(unsigned int i(0); i < size; ++i)
		{
			dTarget[i * nComponents] = d1[i];
		}
	}

	#define BOOST_TT_rep_expression(r, data, T) \
		template void combineArrays(T* d1, unsigned int size, T* dTarget, unsigned int nComponents=1);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
	
	template <typename T> void combineArrays(T* d1,
	                                         T* d2,
	                                         unsigned int size,
	                                         T* dTarget,
	                                         unsigned int nComponents = 2)
	{
		if (nComponents < 2)
			errorMessage("combineArrays() - attempt to provide nComponents that is less than 2");

		for(unsigned int i(0); i < size; ++i)
		{
			dTarget[nComponents * i] = d1[i];
			dTarget[nComponents * i + 1] = d2[i];
		}
	}

	#define BOOST_TT_rep_expression(r, data, T) \
			template void combineArrays(T* d1, T* d2, unsigned int size, T* dTarget, unsigned int nComponents=2);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> void combineArrays(T* d1,
	                                         T* d2,
	                                         T* d3,
	                                         unsigned int size,
	                                         T* dTarget,
	                                         unsigned int nComponents = 3)
	{
		if (nComponents < 3)
			errorMessage("combineArrays() - attempt to provide nComponents that is less than 3");

		for (unsigned int i(0); i < size; ++i)
		{
			dTarget[nComponents * i] = d1[i];
			dTarget[nComponents * i + 1] = d2[i];
			dTarget[nComponents * i + 2] = d3[i];
		}
	}

	#define BOOST_TT_rep_expression(r, data, T) \
		template void combineArrays(T* d1, T* d2, T* d3, unsigned int size, T* dTarget, unsigned int nComponents=3);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
	
	template <typename T> T* combineArrays(T* d1, T* d2, unsigned int size)
	{
		T* d(new T[size * 2]);
		combineArrays(d1, d2, size, d);
		return d;
	}

	#define BOOST_TT_rep_expression(r, data, T) \
		template T* combineArrays(T* d1, T* d2, unsigned int size);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> T* combineArrays(T* d1, T* d2, T* d3, unsigned int size)
	{
		T* d(new T[size * 3]);
		combineArrays(d1, d2, d3, size, d);
		return d;
	}

	#define BOOST_TT_rep_expression(r, data, T) \
		template T* combineArrays(T* d1, T* d2, T* d3, unsigned int size);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename TO, typename TI> TO* combineArrays(TI* d0, TI* d1, TI* d2, TI* d3, unsigned int size)
	{
		TO* d(new TO[size*4]);
		for(unsigned int i(0); i < size; ++i)
		{
			d[4 * i] = d0[i];
			d[4 * i + 1] = d1[i];
			d[4 * i + 2] = d2[i];
			d[4 * i + 3] = d3[i];
		}
		return d;
	}

	template vtkIdType* combineArrays<vtkIdType, unsigned int>(unsigned int* d0,
	                                                           unsigned int* d1,
	                                                           unsigned int* d2,
	                                                           unsigned int* d3,
	                                                           unsigned int size);

	template <typename T> class VTKDataArrayClass{};
	template <> class VTKDataArrayClass<cl_double>{public: typedef vtkDoubleArray T;};
	template <> class VTKDataArrayClass<cl_float>{public: typedef vtkFloatArray T;};
	template <> class VTKDataArrayClass<cl_int>{public: typedef vtkIntArray T;};
	template <> class VTKDataArrayClass<cl_uint>{public: typedef vtkUnsignedIntArray T;};
	template <> class VTKDataArrayClass<cl_long>{public: typedef vtkLongLongArray T;};
	
	template <typename T> AVec<T> castVTKVector(AVec<T> a, T fill=0)
	{
		AVec<T> b(3);
		if (nD(a) == 1)
		{
			b[0] = a[0]; b[1] = fill; b[2] = fill;
		}
		if (nD(a) == 2)
		{
			b[0] = a[1]; b[1] = a[0]; b[2] = fill;
		}
		if (nD(a) == 3)
		{
			b[0] = a[2]; b[1] = a[1]; b[2] = a[0];
		}
		return b;
	}
	
	#define BOOST_TT_rep_expression(r, data, T) \
		template AVec<T> castVTKVector(AVec<T> a, T fill);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	/// Combines 2 arrays into 3 (while 3rd one is made of '0')
	/// and puts them in dTarget
	template <typename T> void combine2ArraysInto3(T* d1,
	                                               T* d2,
	                                               unsigned int size,
	                                               T* dTarget)
	{
		for (unsigned int i(0); i < size; ++i)
		{
			dTarget[3 * i] = d2[i];
			dTarget[3 * i + 1] = d1[i];
			dTarget[3 * i + 2] = 0;
		}
	}

	template <typename T> T* combine2ArraysInto3(T* d1, T* d2, unsigned int size)
	{
		T* d(new T[size * 3]);
		combine2ArraysInto3(d1, d2, size, d);
		return d;
	}
	
	
	template<typename T> vtkSmartPointer<vtkDataArray> castVTKDataArray(T *d,
	                                                                    unsigned int np,
	                                                                    unsigned int save,
	                                                                    const string &name)
	{
		typedef typename VTKDataArrayClass<T>::T DT;
		vtkSmartPointer<DT> vtkArray(vtkSmartPointer<DT>::New());
		vtkArray->SetName(name.c_str());
		vtkArray->SetNumberOfComponents(1);
		typedef decltype(DT::GetDataTypeValueMin()) vtkT;
		if(sizeof(vtkT)!=sizeof(T))
		   errorMessage("castVTKDataArray: type conversion is illegal");
		vtkArray->SetArray((vtkT*) d, np, save);
		return vtkArray;
	}

	
	#define BOOST_TT_rep_expression(r, data, T) \
		template vtkSmartPointer<vtkDataArray> castVTKDataArray<T>(T *d, unsigned int np, unsigned int save, const string &name);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression	

	template<typename T> vtkSmartPointer<vtkDataArray> castVTKDataArray(T *d1, T *d2, unsigned int np, const string &name)
	{
		typedef typename VTKDataArrayClass<T>::T DT;
		vtkSmartPointer<DT> vtkArray(vtkSmartPointer<DT>::New());
		vtkArray->SetName(name.c_str());
		vtkArray->SetNumberOfComponents(2);
		typedef decltype(DT::GetDataTypeValueMin()) vtkT;
		if(sizeof(vtkT)!=sizeof(T))
		   errorMessage("castVTKDataArray: type convetion is illegal");
		vtkArray->SetArray((vtkT*)combineArrays(d1, d2, np), 2 * np, 0);
		return vtkArray;
	}

	#define BOOST_TT_rep_expression(r, data, T) \
			template vtkSmartPointer<vtkDataArray> castVTKDataArray(T *d1, T *d2, unsigned int np, const string &name);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression	

	template<typename T> vtkSmartPointer<vtkDataArray> castVTKDataArray(T *d1, T *d2, T *d3, unsigned int np, const string &name)
	{
		typedef typename VTKDataArrayClass<T>::T DT;
		vtkSmartPointer<DT> vtkArray(vtkSmartPointer<DT>::New());
		vtkArray->SetName(name.c_str());
		vtkArray->SetNumberOfComponents(3);
		typedef decltype(DT::GetDataTypeValueMin()) vtkT;
		if(sizeof(vtkT)!=sizeof(T))
		   errorMessage("castVTKDataArray: type convetion is illegal");
		vtkArray->SetArray((vtkT*)combineArrays(d1, d2, d3, np), 3 * np, 0);
		return vtkArray;
	}

	#define BOOST_TT_rep_expression(r, data, T) \
		template vtkSmartPointer<vtkDataArray> castVTKDataArray(T *d1, T *d2, T *d3, unsigned int np, const string &name);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression	

	template<typename T> vtkSmartPointer<vtkDataArray> castVTKDataArray2in3(T *d1, T *d2, unsigned int np, const string &name)
	{
		typedef typename VTKDataArrayClass<T>::T DT;
		vtkSmartPointer<DT> vtkArray(vtkSmartPointer<DT>::New());
		vtkArray->SetName(name.c_str());
		vtkArray->SetNumberOfComponents(3);
		typedef decltype(DT::GetDataTypeValueMin()) vtkT;
		if(sizeof(vtkT)!=sizeof(T))
		   errorMessage("castVTKDataArray: type convetion is illegal");
		vtkArray->SetArray((vtkT*)combine2ArraysInto3(d1, d2, np), 3 * np, 0);
		return vtkArray;
	}

	#define BOOST_TT_rep_expression(r, data, T) \
		template vtkSmartPointer<vtkDataArray> castVTKDataArray2in3(T *d1, T *d2, unsigned int np, const string &name);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression	

	
	vtkSmartPointer<vtkIdTypeArray> castVTKIdTypeArray(unsigned int *d0, unsigned int *d1, unsigned int *d2, unsigned int *d3, unsigned int np, const std::string &name)	
	{
		vtkSmartPointer<vtkIdTypeArray> vtkArray(vtkSmartPointer<vtkIdTypeArray>::New());
		vtkArray->SetName(name.c_str());
		vtkArray->SetNumberOfComponents(1);
		vtkArray->SetArray(combineArrays<vtkIdType, unsigned int>(d0, d1, d2, d3, np), 4 * np, 0);
		return vtkArray;
	}


	vtkSmartPointer<vtkDataArray> castVTKDataArray(acl::Element source, const string & name)
	{
		if (!isMemBlock(source))
			errorMessage ("castVTKDoubleArray(): provided element is not a MemBlock type");

		if (source->getTypeID() == acl::TYPE_DOUBLE)
		{
			double* d = new double[source->getSize()];
			copy(source, d);
			return castVTKDataArray(d, source->getSize(), 0, name);
		}
		if (source->getTypeID() == acl::TYPE_FLOAT)
		{
			float* d = new float[source->getSize()];
			copy(source, d);
			return castVTKDataArray(d, source->getSize(), 0, name);
		}
		if (source->getTypeID() == acl::TYPE_INT)
		{
			int* d = new int[source->getSize()];
			copy(source, d);
			return castVTKDataArray(d, source->getSize(), 0, name);
		}
		return NULL;
	}	


	vtkSmartPointer<vtkImageData> castVTKData(const Block & b)
	{
		vtkSmartPointer<vtkImageData> image(vtkSmartPointer<vtkImageData>::New());
		image->SetSpacing(b.dx,b.dx,b.dx);
		
		image->SetOrigin(&castVTKVector(b.position)[0]);
		image->SetDimensions(&castVTKVector(b.getSize(),1)[0]);
		return image;
	}


	vtkSmartPointer<vtkImageData> castVTKData(double *d, const Block & b, unsigned int save, const std::string &name)
	{
		vtkSmartPointer<vtkImageData> image(castVTKData(b));
		image->GetPointData()->SetScalars(castVTKDataArray(d,
		                                                   productOfElements(b.getSize()),
		                                                   save,
		                                                   name));
		return image;
	}


	vtkSmartPointer<vtkImageData> castVTKData(double *d1, double *d2, const Block & b, const std::string &name)
	{
		vtkSmartPointer<vtkImageData> image(castVTKData(b));
		image->GetPointData()->SetVectors(castVTKDataArray(d1, d2,
		                                                   productOfElements(b.getSize()),
		                                                   name));
		return image;
	}


	vtkSmartPointer<vtkImageData> castVTKData(double *d1, double *d2, double *d3, const Block & b, const std::string &name)
	{
		vtkSmartPointer<vtkImageData> image(castVTKData(b));
		image->GetPointData()->SetVectors(castVTKDataArray(d1, d2, d3,
		                                                   productOfElements(b.getSize()),
		                                                   name));
		return image;
	}


	vtkSmartPointer<vtkImageData> castVTKData(const AbstractData & d, const vector<string> & names)
	{
		const Block &b(d.getBlock());
		vtkSmartPointer<vtkImageData> image(vtkSmartPointer<vtkImageData>::New());
		image->SetSpacing(b.dx, b.dx, b.dx);
		
		image->SetOrigin(&castVTKVector(b.position)[0]);
		image->SetDimensions(&castVTKVector(b.getSize(), 1)[0]);
		bool bnames(names.size() > 0);
		image->GetPointData()->SetScalars(castVTKDataArray(d.getDContainer()[0], bnames ? names[0] : ""));
		for (unsigned int i(1); i < d.getDContainer().size(); ++i)
			image->GetPointData()->AddArray(castVTKDataArray(d.getDContainer()[i], bnames ? names[i] : ""));
		return image;				
	}


	void putToVTKData(double *d, vtkSmartPointer<vtkImageData> target)
	{
		double *dTarget(((vtkDoubleArray *) target->GetPointData()->GetArray(0))->GetPointer(0));
		unsigned int size(target->GetPointData()->GetArray(0)->GetNumberOfTuples());
		unsigned int nComponents(target->GetPointData()->GetArray(0)->GetNumberOfComponents());
		combineArrays(d, size, dTarget, nComponents);
	}


	void putToVTKData(double *d1, double *d2, vtkSmartPointer<vtkImageData> target)
	{
		double *dTarget(((vtkDoubleArray *) target->GetPointData()->GetArray(0))->GetPointer(0));
		unsigned int size(target->GetPointData()->GetArray(0)->GetNumberOfTuples());
		unsigned int nComponents(target->GetPointData()->GetArray(0)->GetNumberOfComponents());
		combineArrays(d1, d2, size, dTarget, nComponents);
	}


	void putToVTKData(double *d1, double *d2, double *d3, vtkSmartPointer<vtkImageData> target)
	{
		double *dTarget(((vtkDoubleArray *) target->GetPointData()->GetArray(0))->GetPointer(0));
		unsigned int size(target->GetPointData()->GetArray(0)->GetNumberOfTuples());
		unsigned int nComponents(target->GetPointData()->GetArray(0)->GetNumberOfComponents());
		combineArrays(d1, d2, d3, size, dTarget, nComponents);
	}


	template <typename T> void decomposeArrays(T* d, T* d1, T* d2, unsigned int size)
	{
		for (unsigned int i(0); i < size; ++i)
		{
			d1[i] = d[2 * i];
			d2[i] = d[2 * i + 1];
		}
	}

	template void decomposeArrays<double>(double* d, double* d1, double* d2, unsigned int size);
	template void decomposeArrays<float>(float* d, float* d1, float* d2, unsigned int size);
	template void decomposeArrays<int>(int* d, int* d1, int* d2, unsigned int size);
	template void decomposeArrays<unsigned int>(unsigned int* d, unsigned int* d1, unsigned int* d2, unsigned int size);

	template <typename T> void decomposeArrays(T* d, T* d1, T* d2, T* d3, unsigned int size)
	{
		for (unsigned int i(0); i < size; ++i)
		{
			d1[i] = d[3 * i];
			d2[i] = d[3 * i + 1];
			d3[i] = d[3 * i + 2];
		}
	}

	template void decomposeArrays<double>(double* d, double* d1, double* d2, double* d3, unsigned int size);
	template void decomposeArrays<float>(float* d, float* d1, float* d2, float* d3, unsigned int size);
	template void decomposeArrays<int>(int* d, int* d1, int* d2, int* d3, unsigned int size);
	template void decomposeArrays<unsigned int>(unsigned int* d, unsigned int* d1, unsigned int* d2, unsigned int* d3, unsigned int size);


	// Creates destination arrays (they need to be deleted outside after being used!),
	// transfers data from \p source to \p destination.
	template <typename T> void decomposeVTKDataArray(vtkDataArray * source,
	                                                 vector<T *> & destination)
	{
		typedef typename VTKDataArrayClass<T>::T DataType;
		
		unsigned int numOfComponents(source->GetNumberOfComponents());
		unsigned int numOfTuples(source->GetNumberOfTuples());

		for (unsigned int i = 0; i < numOfComponents; ++i)
		{
			destination.push_back(new T[numOfTuples]);
		}

		DataType * array((DataType *) source);
		for (unsigned int i = 0; i < numOfTuples; ++i)
		{
			for (unsigned int j = 0; j < numOfComponents; ++j)
			{
				destination[j][i] = array->GetValue(numOfComponents * i + j);
			}
		}
	}

	template void decomposeVTKDataArray<double>(vtkDataArray * source,
	                                            vector<double *> & destination);
	template void decomposeVTKDataArray<float>(vtkDataArray * source,
	                                           vector<float *> & destination);
	template void decomposeVTKDataArray<int>(vtkDataArray * source,
	                                         vector<int *> & destination);
	template void decomposeVTKDataArray<unsigned int>(vtkDataArray * source,
	                                                  vector<unsigned int *> & destination);

	std::shared_ptr<Block> makeBlock(vtkSmartPointer<vtkImageData> image)
	{
		int dims[3];
		image->GetDimensions(dims);

		double spacing[3];
		image->GetSpacing(spacing);

		double origin[3];
		int extent[6];
		image->GetOrigin(origin);
		image->GetExtent(extent);
		origin[0]+=extent[0]*spacing[0];
		origin[1]+=extent[2]*spacing[1];
		origin[2]+=extent[4]*spacing[2];
		
		shared_ptr<Block> block(new Block(makeAVec<int>(dims[2], dims[1], dims[0]),
		                                  spacing[0],
		                                  makeAVec<double>(origin[2], origin[1], origin[0])));

		return block;
	}


	// Converts vtk data type into asl data type;
	// exits if non-supported data type is provided.
	acl::TypeID aslType(int vtkType)
	{
		acl::TypeID t(acl::TYPE_INT);
		switch (vtkType)
		{
			case VTK_INT :
				t = acl::TYPE_INT;
				break;
			case VTK_UNSIGNED_INT :
				t = acl::TYPE_UINT;
				break;
			case VTK_FLOAT :
				t = acl::TYPE_FLOAT;
				break;
			case VTK_DOUBLE :
				t = acl::TYPE_DOUBLE;
				break;
			case VTK_LONG :
				t = acl::TYPE_LONG;
				break;
			default :
				errorMessage("aslType - vtk data type not recognized: " + numToStr(vtkType));
				break;
		}
		return t;
	}
	
	SPDataWithGhostNodesACLData makeData(vtkSmartPointer<vtkImageData> image,
	                                     unsigned int arrayNum,
	                                     acl::CommandQueue queue)
	{

		if ((int)arrayNum >= image->GetPointData()->GetNumberOfArrays())
			errorMessage("makeData() - arrayNum out of range");

		shared_ptr<Block> block = makeBlock(image);
		acl::TypeID type = aslType(image->GetPointData()->GetArray(arrayNum)->GetDataType());
		unsigned int numOfComponents(image->GetPointData()->GetArray(arrayNum)->GetNumberOfComponents());
		SPDataWithGhostNodesACLData data = generateDataContainerACL_SP(offset(*block,-1),
		                                                               type,
		                                                               numOfComponents,
		                                                               1,
		                                                               queue);

		switch (type)
		{
			case acl::TYPE_INT :
			{
				vector<cl_int *> decomposedData;
				decomposeVTKDataArray<cl_int>(image->GetPointData()->GetArray(arrayNum),
				                              decomposedData);
				for (unsigned int i = 0; i < decomposedData.size(); ++i)
				{
					copy(decomposedData[i], data->getDContainer()[i]);
					delete[] decomposedData[i];
				}
				break;
			}
			case acl::TYPE_UINT :
			{
				vector<cl_uint *> decomposedData;
				decomposeVTKDataArray<cl_uint>(image->GetPointData()->GetArray(arrayNum),
				                               decomposedData);
				for (unsigned int i = 0; i < decomposedData.size(); ++i)
				{
					copy(decomposedData[i], data->getDContainer()[i]);
					delete[] decomposedData[i];
				}
				break;
			}
			case acl::TYPE_FLOAT :
			{
				vector<cl_float *> decomposedData;
				decomposeVTKDataArray<cl_float>(image->GetPointData()->GetArray(arrayNum),
				                                decomposedData);
				for (unsigned int i = 0; i < decomposedData.size(); ++i)
				{
					copy(decomposedData[i], data->getDContainer()[i]);
					delete[] decomposedData[i];
				}
				break;
			}
			case acl::TYPE_DOUBLE :
			{
				vector<cl_double *> decomposedData;
				decomposeVTKDataArray<cl_double>(image->GetPointData()->GetArray(arrayNum),
				                                 decomposedData);
				for (unsigned int i = 0; i < decomposedData.size(); ++i)
				{
					copy(decomposedData[i], data->getDContainer()[i]);
					delete[] decomposedData[i];
				}
				break;
			}
			case acl::TYPE_LONG :
			{
				vector<cl_long *> decomposedData;
				decomposeVTKDataArray<cl_long>(image->GetPointData()->GetArray(arrayNum),
				                               decomposedData);
				for (unsigned int i = 0; i < decomposedData.size(); ++i)
				{
					copy(decomposedData[i], data->getDContainer()[i]);
					delete[] decomposedData[i];
				}
				break;
			}
		}		

		return data;
	}
	
		
} // asl

