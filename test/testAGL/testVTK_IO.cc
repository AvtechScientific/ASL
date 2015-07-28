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


/**
	\example testVTK_IO.cc
 */

#include <writers/aslVTKFormatWriters.h>
#include <readers/aslVTKFormatReaders.h>
#include <aslGenerators.h>
#include <num/aslDataResampling.h>
#include <data/aslDataWithGhostNodes.h>

void testMINC()
{
	cout << "Test of MINC files reader..." << flush;

	auto data(asl::read("subject04_crisp_v.mnc", 0));
	
	asl::writeVTKXML("data.vti",
	                 *data,
	                 "data");
}

void testMINCplus()
{
	cout << "Test of MINC files reader +..." << flush;

	auto data(asl::read("subject04_crisp_v.mnc", 0));
	
	asl::DataCoarser dc(data);
	dc.init();
	dc.execute();
	asl::writeVTKXML("dataCoarsed.vti",
	                 *dc.getDataOut(),
	                 "data");

}

void testSurfaceSTL()
{
	cout << "Test of Surface STL files reader..." << flush;

	auto data(asl::readSurface("bus.stl", 5));
//	auto data(asl::readSurface("xx.vtp", .01));
	
	asl::writeVTKXML("dataSurfaceSTL.vti",
	                 *data,
	                 "data");
}



int main()
{
//	testMINC();
//	testMINCplus();
	testSurfaceSTL();
	return 0;
}
