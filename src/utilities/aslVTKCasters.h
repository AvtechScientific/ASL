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


#ifndef ASLVTKCASTERS_H
#define ASLVTKCASTERS_H

#include <memory>
#include <string>
#include <vector>
#include <vtkSmartPointer.h>
#include <acl/aclHardware.h>
#include <aslGenerators.h>
#include <math/aslVectors.h>

/**
	\defgroup VTKInterfacing VTK Interfacing
	\ingroup Interfacing
*/

class vtkDataArray; 
class vtkImageData;
class vtkIdTypeArray;

namespace acl
{
	class ElementBase;
	typedef std::shared_ptr<ElementBase> Element;
}

namespace asl
{
	class Block;
	class AbstractData;
	
	
	/// @{
	/// \ingroup VTKInterfacing

	/// creates VTKDataArray with 1 component \p d and length \p np and \p name 
	template<typename T> vtkSmartPointer<vtkDataArray> castVTKDataArray(T *d,
	                                                                    unsigned int np,
	                                                                    unsigned int save = 0,
	                                                                    const std::string &name = "");

	/// creates VTKDataArray with 2 component \p d1 and \p d2 and length \p np and \p name 
	template<typename T> vtkSmartPointer<vtkDataArray> castVTKDataArray(T *d1,
	                                                                    T *d2,
	                                                                    unsigned int np,
	                                                                    const std::string &name = "");	

	/// creates VTKDataArray with 3 component \p d1, \p d2 and \p d3 and length \p np and \p name 
	template<typename T> vtkSmartPointer<vtkDataArray> castVTKDataArray(T *d1,
	                                                                    T *d2,
	                                                                    T *d3,
	                                                                    unsigned int np,
	                                                                    const std::string &name = "");	

	/// creates VTKDataArray with 3 component \p d2, \p d1 and 0 and length \p np and \p name 
	template<typename T> vtkSmartPointer<vtkDataArray> castVTKDataArray2in3(T *d1, T *d2, 
	                                                                        unsigned int np, 
	                                                                        const std::string &name);

	
	/// creates VTKDataArray with 3 component \p d1, \p d2 and \p d3 and length \p np and \p name 
	vtkSmartPointer<vtkIdTypeArray> castVTKIdTypeArray(unsigned int *d0,
	                                                   unsigned int *d1,
	                                                   unsigned int *d2,
	                                                   unsigned int *d3,
	                                                   unsigned int np,
	                                                   const std::string &name = "");	

	vtkSmartPointer<vtkDataArray> castVTKDataArray(acl::Element source,
	                                               const std::string &name = "");

	vtkSmartPointer<vtkImageData> castVTKData(const Block & b);	

	vtkSmartPointer<vtkImageData> castVTKData(double *d,
	                                          const Block & b,
	                                          unsigned int save = 0,
	                                          const std::string &name = "");

	vtkSmartPointer<vtkImageData> castVTKData(double *d1,
	                                          double *d2,
	                                          const Block & b,
	                                          const std::string &name = "");

	vtkSmartPointer<vtkImageData> castVTKData(double *d1,
	                                          double *d2,
	                                          double *d3,
	                                          const Block & b,
	                                          const std::string &name = "");

	void putToVTKData(double *d, vtkSmartPointer<vtkImageData> target);
	void putToVTKData(double *d1, double *d2, vtkSmartPointer<vtkImageData> target);
	void putToVTKData(double *d1, double *d2, double *d3, vtkSmartPointer<vtkImageData> target);

	
	vtkSmartPointer<vtkImageData> castVTKData(const AbstractData & d,
	                                          const std::vector<std::string> &names = std::vector<std::string>(0));	

	SPDataWithGhostNodesACLData makeData(vtkSmartPointer<vtkImageData> image,
	                                     unsigned int arrayNum = 0,
	                                     acl::CommandQueue queue = acl::hardware.defaultQueue);
	
	std::shared_ptr<Block> makeBlock(vtkSmartPointer<vtkImageData> image);
	
	template<typename T> inline vtkSmartPointer<vtkDataArray> castVTKDataArray(std::shared_ptr<T> d,
	                                                                           unsigned int np,
	                                                                           const std::string &name = "")
	{
		return castVTKDataArray(d.get(), np, 1, name);	
	}
	
	vtkSmartPointer<vtkImageData> inline castVTKData(std::shared_ptr<double> d,
	                                                 const Block & b,
	                                                 const std::string &name = "")
	{
		return castVTKData(d.get(), b, 1, name);	
	}


	template<typename T> inline vtkSmartPointer<vtkDataArray> castVTKDataArray(std::vector<T> & d,
	                                                                           unsigned int np,
	                                                                           const std::string &name = "")
	{
		return castVTKDataArray(&d.front(), np, 1, name);	
	}
	

	vtkSmartPointer<vtkImageData> inline castVTKData(std::vector<double> & d,
	                                                 const Block & b,
	                                                 const std::string &name = "")
	{
		return castVTKData(&d.front(), b, 1, name);	
	}


	template <typename T> AVec<T> castVTKVector(AVec<T> a, T fill = 0);


	template <typename T> void combineArrays(T* d1,
	                                         T* d2,
	                                         T* d3,
	                                         unsigned int size,
	                                         T* dTarget,
	                                         unsigned int nComponents = 3);
	/// @}
	
}  //namespace acl

//--------------------------------Implementation------------------------


#endif // ASLVTKCASTERS_H
