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


#include "aslVTKFormatWriters.h"
#include <utilities/aslVTKCasters.h>
//#include <agl/vtk/aglVTK.h>

#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkImageData.h>

#include "math/aslVectors.h"
#include "acl/acl.h"
#include <acl/DataTypes/aclMemBlock.h>
#include <acl/aclTypesList.h>
#include <data/aslDataWrapper.h>

namespace asl
{

	vtkSmartPointer<vtkImageData> makeVTKData(const Block & block,
	                                          const vector<pair<string, acl::VectorOfElementsData> > & scalarFields,
	                                          const vector<pair<string, acl::VectorOfElementsData> > & vectorFields)
	{
		vtkSmartPointer<vtkImageData> image(vtkSmartPointer<vtkImageData>::New());
		image->SetSpacing(block.dx, block.dx, block.dx);
		
		image->SetOrigin(&castVTKVector(block.position)[0]);
		image->SetDimensions(&castVTKVector(block.getSize(), 1)[0]);

		if ((scalarFields.size() == 0) && (vectorFields.size() == 0))
			errorMessage("WriterVTK::makeVTKData() - there are no fields to write");

		for (unsigned int i = 0; i < scalarFields.size(); ++i)
		{
			auto vtkArray(castVTKDataArray(scalarFields[i].second[0],
			                                scalarFields[i].first));
			if (i==0)
				image->GetPointData()->SetScalars(vtkArray);
			image->GetPointData()->AddArray(vtkArray);
		}

		for (unsigned int i = 0; i < vectorFields.size(); ++i)
		{
			// compose a long vtk vector array and then feed it below
			vtkSmartPointer<vtkDataArray> vtkArray;
			switch (vectorFields[i].second[0]->getTypeID()) 
			{
				#define BOOST_TT_rep_expression(r, data, t) \
					case acl::typeToTypeID<t>(): \
					{ \
						if (vectorFields[i].second.size() == 2) \
						{   \
							auto p0(acl::map<t>(vectorFields[i].second[0])); \
							auto p1(acl::map<t>(vectorFields[i].second[1])); \
							vtkArray = castVTKDataArray2in3(p0.get(), \
												            p1.get(), \
												            vectorFields[i].second[0]->getSize(), \
												            vectorFields[i].first); \
						}   \
						if (vectorFields[i].second.size() == 3) \
						{   \
							auto p0(acl::map<t>(vectorFields[i].second[2])); \
							auto p1(acl::map<t>(vectorFields[i].second[1])); \
							auto p2(acl::map<t>(vectorFields[i].second[0])); \
							vtkArray = castVTKDataArray(p0.get(), \
												        p1.get(), \
												        p2.get(), \
												        vectorFields[i].second[0]->getSize(), \
												        vectorFields[i].first); \
						}   \
					} \
					break;	
				BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
				#undef BOOST_TT_rep_expression	
			}
			
			image->GetPointData()->AddArray(vtkArray); 
		}	
	
		return image;
	}


	WriterVTKXML::WriterVTKXML(const string & file, Block *nbl) :
		Writer(file),
		newBl(nbl)
	{
	}

	void WriterVTKXML::write()
	{
		vtkSmartPointer<vtkXMLImageDataWriter> writer(vtkSmartPointer<vtkXMLImageDataWriter>::New());
		writer->SetInputData(makeVTKData(*block, scalarFields, vectorFields).GetPointer());
		writer->SetFileName((file + "_" + numToStr(numOfWrites) + ".vti").c_str());
		writer->SetDataModeToAppended();
		writer->EncodeAppendedDataOff();		
//		writer->SetDataModeToBinary();
//		writer->SetDataModeToAscii();
		writer->Write();

		++numOfWrites;
	}


	void writeVTKXML(const string & fileName,
	                 const AbstractData & data,
	                 const string & name)
	{
		vtkSmartPointer<vtkXMLImageDataWriter> writer(vtkSmartPointer<vtkXMLImageDataWriter>::New());
		unsigned int nComp(data.getDContainer().size());

		std::vector<string> names(nComp, name);
		for (unsigned int i(0); i < nComp; ++i)
			names[i] += "-" + numToStr(i);
		
		writer->SetInputData(castVTKData(data, names).GetPointer());
		writer->SetFileName(fileName.c_str());
		writer->SetDataModeToAppended();
		writer->EncodeAppendedDataOff();		
//		writer->SetDataModeToBinary();
//		writer->SetDataModeToAscii();
		writer->Write();
	}
	
} //asl


