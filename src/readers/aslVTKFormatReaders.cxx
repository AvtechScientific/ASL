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


#include "aslVTKFormatReaders.h"
#include <acl/acl.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <utilities/aslVTKCasters.h>
#include <vtkMINCImageReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkImplicitModeller.h>
#include <boost/filesystem.hpp>
#include <data/aslDataWithGhostNodes.h>

using namespace boost::filesystem;

namespace asl
{
	using namespace asl;
	
	SPDataWithGhostNodesACLData readMINC(const string & fileName,
	                                     unsigned int arrayNum,
	                                     acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkMINCImageReader> reader(vtkSmartPointer<vtkMINCImageReader>::New());
		reader->RescaleRealValuesOn();
		if (!reader->CanReadFile(fileName.c_str()))
			errorMessage("MINC reader: The input file is corrupted or file name is wrong: " + fileName);
		reader->SetFileName(fileName.c_str());
		reader->Update();
		auto data(makeData(reader->GetOutput(), arrayNum, queue));
		return data;
	}


	SPDataWithGhostNodesACLData readVTKXML(const string & fileName,
	                                       unsigned int arrayNum,
	                                       acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkXMLImageDataReader> reader(vtkSmartPointer<vtkXMLImageDataReader>::New());
		if (!reader->CanReadFile(fileName.c_str()))
			errorMessage("VTKXML reader: The input file is corrupted or file name is wrong: " + fileName);
		reader->SetFileName(fileName.c_str());
		reader->Update();
		auto data(makeData(reader->GetOutput(), arrayNum, queue));
		return data;
	}


	SPDataWithGhostNodesACLData readVTK(const string & fileName,
	                                    unsigned int arrayNum,
	                                    acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkStructuredPointsReader> reader(vtkSmartPointer<vtkStructuredPointsReader>::New());
		reader->SetFileName(fileName.c_str());
		reader->Update();
		if (!reader->IsFileStructuredPoints())
			errorMessage("VTK reader: The input file is corrupted or file name is wrong: " + fileName);

		auto data(makeData(vtkSmartPointer<vtkStructuredPoints>(reader->GetOutput()), arrayNum, queue));
		return data;
	}


	SPDataWithGhostNodesACLData readDICOM(const string & fileName,
	                                      unsigned int arrayNum,
	                                      acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkDICOMImageReader> reader(vtkSmartPointer<vtkDICOMImageReader>::New());
		if (!reader->CanReadFile(fileName.c_str()))
			errorMessage("DICOM reader: The input file is corrupted or file name is wrong: " + fileName);
		reader->SetFileName(fileName.c_str());
		reader->Update();
		auto data(makeData(reader->GetOutput(), arrayNum, queue));
		return data;
	}


	SPDataWithGhostNodesACLData read(const string & fileName,
	                                 unsigned int arrayNum,
	                                 acl::CommandQueue queue)
	{
		string fileExtension;
		path pathToFile(fileName);
		fileExtension = pathToFile.extension().string();

		SPDataWithGhostNodesACLData data;

		if (fileExtension == ".vtk")
			data = readVTK(fileName, arrayNum, queue);

		if (fileExtension == ".vti")
			data = readVTKXML(fileName, arrayNum, queue);

		if (fileExtension == ".mnc")
			data = readMINC(fileName, arrayNum, queue);

		if (fileExtension == ".dcm")
			data = readDICOM(fileName, arrayNum, queue);

		if (data.get() == 0)
			errorMessage("Reader: file format not supported");

		return data;
	}

	SPDataWithGhostNodesACLData surfaceToData(vtkDataSet* surf,
	                                          double dx,
	                                          acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkImplicitModeller> modeller(vtkSmartPointer<vtkImplicitModeller>::New());
		modeller->SetOutputScalarTypeToFloat();
		modeller->CappingOff();
		modeller->SetCapValue(dx*3.6);
		modeller->SetInputData(surf);
		modeller->ComputeModelBounds(surf);
		modeller->SetMaximumDistance(dx*1.8*2.);
		double b[6];
		modeller->GetModelBounds (b);
		modeller->SetSampleDimensions((b[1]-b[0]) / dx, (b[3]-b[2]) / dx, (b[5]-b[4]) / dx);
		modeller->SetProcessModeToPerVoxel();
		modeller->Update();

		auto data(makeData(modeller->GetOutput(), 0, queue));
		acl::initData(data->getEContainer(), (data->getEContainer() - dx*1.8)/dx/1.8);

		return data;
	}

	SPDataWithGhostNodesACLData surfaceToData(vtkDataSet* surf,
	                                          Block & b,
	                                          acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkImplicitModeller> modeller(vtkSmartPointer<vtkImplicitModeller>::New());
		modeller->SetOutputScalarTypeToFloat();
		modeller->CappingOff();
		modeller->SetCapValue(b.dx*3.6);
		modeller->SetInputData(surf);
		modeller->SetMaximumDistance(b.dx*1.8*2.);
		double bound[6];
		auto bp(b.getBPosition()); 
		bound[0]=b.position[2]; bound[1]=bp[2];
		bound[2]=b.position[1]; bound[3]=bp[1];
		bound[4]=b.position[0]; bound[5]=bp[0];		
		modeller->SetModelBounds (bound);
		int s[3];
		s[0] = b.getSize()[2];
		s[1] = b.getSize()[1];
		s[2] = b.getSize()[0];
		modeller->SetSampleDimensions(s);
		modeller->SetProcessModeToPerVoxel();
		modeller->Update();

		auto data(makeData(modeller->GetOutput(), 0, queue));
		acl::initData(data->getEContainer(), (data->getEContainer() - b.dx*1.2)/b.dx/1.8);

		return data;
	}

	SPDataWithGhostNodesACLData surfaceToData(vtkDataSet* surf,
	                                          double dx,
	                                          double offset_X0, double offset_XE,
	                                          double offset_Y0, double offset_YE,
	                                          double offset_Z0, double offset_ZE,
	                                          acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkImplicitModeller> modeller(vtkSmartPointer<vtkImplicitModeller>::New());
		modeller->SetOutputScalarTypeToFloat();
		modeller->CappingOff();
		modeller->SetCapValue(dx*3.6);
		modeller->SetInputData(surf);
		modeller->ComputeModelBounds(surf);
		modeller->SetMaximumDistance(dx*1.8*2.);
		double b[6];
		modeller->GetModelBounds (b);
		double bound[6];		
		bound[0]=b[0] - (b[1]-b[0])*offset_Z0; bound[1]=b[1] + (b[1]-b[0])*offset_ZE;
		bound[2]=b[2] - (b[3]-b[2])*offset_Y0; bound[3]=b[3] + (b[3]-b[2])*offset_YE;
		bound[4]=b[4] - (b[5]-b[4])*offset_X0; bound[5]=b[5] + (b[5]-b[4])*offset_XE;
		modeller->SetModelBounds (bound);
		modeller->SetSampleDimensions((bound[1]-bound[0]) / dx, 
		                              (bound[3]-bound[2]) / dx, 
		                              (bound[5]-bound[4]) / dx);
		modeller->SetProcessModeToPerVoxel();
		modeller->Update();

		auto data(makeData(modeller->GetOutput(), 0, queue));
		acl::initData(data->getEContainer(), (data->getEContainer() - dx*1.8)/dx/1.8);

		return data;
	}

	
	SPDataWithGhostNodesACLData readSurfaceVTKXML(const string & fileName,
	                                              double dx,
	                                              acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkXMLPolyDataReader> reader(vtkSmartPointer<vtkXMLPolyDataReader>::New());
		if (!reader->CanReadFile(fileName.c_str()))
			errorMessage("SurfaceVTKXML reader: The input file is corrupted or file name is wrong: " + fileName);
		reader->SetFileName(fileName.c_str());
		reader->Update();

		return surfaceToData(reader->GetOutput(),dx,queue);
	}

	SPDataWithGhostNodesACLData readSurfaceVTKXML(const string & fileName,
	                                              Block & b,
	                                              acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkXMLPolyDataReader> reader(vtkSmartPointer<vtkXMLPolyDataReader>::New());
		if (!reader->CanReadFile(fileName.c_str()))
			errorMessage("SurfaceVTKXML reader: The input file is corrupted or file name is wrong: " + fileName);
		reader->SetFileName(fileName.c_str());
		reader->Update();

		return surfaceToData(reader->GetOutput(),b,queue);
	}

	SPDataWithGhostNodesACLData readSurfaceVTKXML(const string & fileName,
	                                              double dx,
	                                              double offset_X0, double offset_XE,
	                                              double offset_Y0, double offset_YE,
	                                              double offset_Z0, double offset_ZE,
	                                              acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkXMLPolyDataReader> reader(vtkSmartPointer<vtkXMLPolyDataReader>::New());
		if (!reader->CanReadFile(fileName.c_str()))
			errorMessage("SurfaceVTKXML reader: The input file is corrupted or file name is wrong: " + fileName);
		reader->SetFileName(fileName.c_str());
		reader->Update();

		return surfaceToData(reader->GetOutput(),
		                     dx,
		                     offset_X0, offset_XE,
		                     offset_Y0, offset_YE,
		                     offset_Z0, offset_ZE,
		                     queue);
	}
	
	SPDataWithGhostNodesACLData readSurfaceSTL(const string & fileName,
	                                           double dx,
	                                           acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkSTLReader> reader(vtkSmartPointer<vtkSTLReader>::New());
		reader->SetFileName(fileName.c_str());
		reader->Update();
		
		return surfaceToData(reader->GetOutput(),dx,queue);
	}
	
	SPDataWithGhostNodesACLData readSurfaceSTL(const string & fileName,
	                                           Block & b,
	                                           acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkSTLReader> reader(vtkSmartPointer<vtkSTLReader>::New());
		reader->SetFileName(fileName.c_str());
		reader->Update();
		
		return surfaceToData(reader->GetOutput(),b,queue);
	}

	SPDataWithGhostNodesACLData readSurfaceSTL(const string & fileName,
	                                           double dx,
	                                           double offset_X0, double offset_XE,
	                                           double offset_Y0, double offset_YE,
	                                           double offset_Z0, double offset_ZE,
	                                           acl::CommandQueue queue)
	{
		vtkSmartPointer<vtkSTLReader> reader(vtkSmartPointer<vtkSTLReader>::New());
		reader->SetFileName(fileName.c_str());
		reader->Update();
		
		return surfaceToData(reader->GetOutput(),
		                     dx,
		                     offset_X0, offset_XE,
		                     offset_Y0, offset_YE,
		                     offset_Z0, offset_ZE,
		                     queue);
	}
	
	SPDataWithGhostNodesACLData readSurface(const string & fileName,
	                                        double dx,
	                                        acl::CommandQueue queue)
	{
		string fileExtension;
		path pathToFile(fileName);
		fileExtension = pathToFile.extension().string();

		SPDataWithGhostNodesACLData data;

		if (fileExtension == ".vtp")
			data = readSurfaceVTKXML(fileName, dx, queue);

		if (fileExtension == ".stl")
			data = readSurfaceSTL(fileName, dx, queue);
		
		if (data.get() == 0)
			errorMessage("Reader: file format not supported");

		return data;		
	}
	
	SPDataWithGhostNodesACLData readSurface(const string & fileName,
	                                        Block & b,
	                                        acl::CommandQueue queue)
	{
		string fileExtension;
		path pathToFile(fileName);
		fileExtension = pathToFile.extension().string();

		SPDataWithGhostNodesACLData data;

		if (fileExtension == ".vtp")
			data = readSurfaceVTKXML(fileName, b, queue);

		if (fileExtension == ".stl")
			data = readSurfaceSTL(fileName, b, queue);
		
		if (data.get() == 0)
			errorMessage("Reader: file format not supported");

		return data;		
	}

	SPDataWithGhostNodesACLData readSurface(const string & fileName,
	                                        double dx,
	                                        double offset_X0, double offset_XE,
	                                        double offset_Y0, double offset_YE,
	                                        double offset_Z0, double offset_ZE,
	                                        acl::CommandQueue queue)
	{
		string fileExtension;
		path pathToFile(fileName);
		fileExtension = pathToFile.extension().string();

		SPDataWithGhostNodesACLData data;

		if (fileExtension == ".vtp")
			data = readSurfaceVTKXML(fileName,
			                         dx,
			                         offset_X0, offset_XE,
			                         offset_Y0, offset_YE,
			                         offset_Z0, offset_ZE,
			                         queue);

		if (fileExtension == ".stl")
			data = readSurfaceSTL(fileName,
			                      dx,
			                      offset_X0, offset_XE,
			                      offset_Y0, offset_YE,
			                      offset_Z0, offset_ZE,
			                      queue);
		
		if (data.get() == 0)
			errorMessage("Reader: file format not supported");

		return data;		
	}
	
} //asl


