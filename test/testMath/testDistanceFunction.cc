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
	\example testDistanceFunction.cc
 */

#include "math/aslVectors.h"
#include "aslGenerators.h"
#include "writers/aslVTKFormatWriters.h"
#include "aslGeomInc.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include "data/aslDataWithGhostNodes.h"
#include "num/aslDFOptimizer.h"
#include "math/aslTemplates.h"

//typedef float FlT;
typedef double FlT;

using asl::AVec;
using asl::makeAVec;

bool testDistFOperations2D()
{
	// Geometry description
	// Radius
	FlT r(10.);

	// Generates a sphere with radius r and center at (50., 50.)
	auto df1(generateDFSphere(r, asl::makeAVec(50., 50.)));
	auto df2(generateDFSphere(r, asl::makeAVec(40., 40.)));
	auto df3(generateDFSphere(2. * r, asl::makeAVec(50., 50.)));
	// Resulting geometry: union of the spheres df1 and df2 intersected by df3 
	auto resultGeometry((df1 | df2) & df3);

	// Geometry to Data conversion
	// Grid size (= discrete size of the simulated domain) 
	asl::AVec<int> size(asl::makeAVec(100., 100.));
	// Grid resolution (= space step)
	FlT dx(1.);
	// Creates a Block which describes the grid
	asl::Block block(size, dx);
	// Allocates memory for the data that corresponds to
	// the nodes of the grid desribed by the \p block.
	auto data(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	// Initializes the \p data with the values of the distance function
	// extracted from all points of the \p resultGeometry.
	asl::initData(data, resultGeometry);

	// Writes the \p data into the file.
	asl::writeVTKXML("distFOperation2D.vti", *data, "data");

	return true;
}


bool testDistFOperations3D()
{
	FlT r(10.);
	FlT dx(1.);
	asl::AVec<int> size(asl::makeAVec(50.,50.,50.));

	asl::Block block(size, dx);
	auto data(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	auto df1(generateDFSphere(r, asl::AVec<FlT>(size)*.5));
	auto df2(generateDFSphere(r, asl::AVec<FlT>(size)*.4));
	auto df3(generateDFSphere(1.5*r, asl::AVec<FlT>(size)*.5));
	asl::initData(data, ((df1 | df2) & df3));

	asl::writeVTKXML("distFOperation3D.vti", *data, "data");

	return true;	
}


bool testDistFOrderedCylinders()
{
	FlT r(3.);
	FlT spacing(4.);
	FlT dx(1.);
	asl::AVec<int> size(asl::makeAVec(50., 50., 50.));

	asl::Block block(size, dx);
	auto data(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	vector<asl::SPDistanceFunction> cylinders;
	asl::SPDistanceFunction resultGeometry;
	asl::AVec<FlT> orientation(asl::makeAVec(0., 0., 1.));
	for (int i = 0; i < size[0] / (2 * r + spacing); ++i)
	{
		for (int j = 0; j < size[1] / (2 * r + spacing); ++j)
		{
			cylinders.push_back(generateDFCylinderInf(r, orientation, asl::makeAVec(i * (2. * r + spacing) + r + spacing / 2., j * (2. * r + spacing) + r + spacing / 2., 0.)));
			resultGeometry = resultGeometry | cylinders.back();
		}
	}

	asl::initData(data, resultGeometry);

	asl::writeVTKXML("distFOrderedCylinders.vti", *data, "data");

	return true;	
}


bool testDistFUnorderedCylinders()
{
	FlT r(3.);
	FlT spacing(4.);
	FlT dx(1.);
	asl::AVec<int> size(asl::makeAVec(100., 100., 100.));

	asl::Block block(size, dx);
	auto data(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	vector<asl::SPDistanceFunction> cylinders;
	asl::SPDistanceFunction resultGeometry;
	asl::AVec<FlT> orientation(asl::makeAVec(0., 0., 1.));
	srand (time(NULL));
	for (unsigned int i = 0; i < size[0] / (2 * r + spacing); ++i)
	{
		for (unsigned int j = 0; j < size[1] / (2 * r + spacing); ++j)
		{
			for (unsigned int d = 0; d < orientation.getSize(); ++d)
				orientation[d] = rand() % size[d];

			cylinders.push_back(generateDFCylinderInf(r, orientation, asl::makeAVec(i * (2. * r + spacing) + r + spacing / 2., j * (2. * r + spacing) + r + spacing / 2., (FlT) (rand() % size[2]))));
			resultGeometry = resultGeometry | cylinders.back();
		}
	}

	asl::initData(data, resultGeometry);

	asl::writeVTKXML("distFUnorderedCylinders.vti", *data, "data");

	return true;	
}


bool testDistFNormalization2D()
{
	FlT r(10.);
	FlT dx(1.);
	asl::AVec<int> size(asl::makeAVec(100.,100.));

	asl::Block block(size, dx);
	auto data(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	auto df1(generateDFSphere(r, asl::AVec<FlT>(size)*.5));
	asl::initData(data, normalize(df1, dx));

	asl::writeVTKXML("distFNormalization2D.vti", *data, "data");

	return true;	
}


bool testDistFNormalization3D()
{
	FlT r(10.);
	FlT dx(1.);
	asl::AVec<int> size(asl::makeAVec(50.,50.,50.));

	asl::Block block(size, dx);
	auto data(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	auto df1(generateDFSphere(r, asl::AVec<FlT>(size)*.5));
	asl::initData(data, normalize(df1,dx));

	asl::writeVTKXML("distFNormalization3D.vti", *data, "data");

	return true;	
}


bool testDistFOperations3DPrism()
{
	FlT r(10.);
	FlT dx(1.);
	asl::AVec<int> size(asl::makeAVec(50.,50.,50.));

	asl::Block block(size, dx);
	auto data(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	auto center(asl::AVec<FlT>(size)*.5);
	auto df1(generateDFSphere(r, center));
	auto df2(asl::generateDFConvexPolygonPrism({center+asl::makeAVec(4.,0.,0.),
                                                              center+asl::makeAVec(4.,4.,0.),
                                                              center+asl::makeAVec(-4.,0.,0.),
                                                              center+asl::makeAVec(-4.,-4.,0.)}));
	asl::initData(data, (df1 & (-df2)));

	asl::writeVTKXML("distFOperation3DPrism.vti", *data, "data");

	return true;	
}


bool testDistFOperations3DBlock()
{
	FlT r(10.);
	FlT dx(1.);
	asl::AVec<int> size(asl::makeAVec(50.,50.,50.));

	asl::Block block(size, dx);
	auto data(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	auto center(asl::AVec<FlT>(size)*.5);
	auto df1(generateDFSphere(r, center));
	auto df2(generateDFInBlock(block, 1));

	asl::initData(data, (df1 | df2));

	asl::writeVTKXML("distFOperation3DBlock.vti", *data, "data");

	return true;	
}

bool testDistFAdvanced3D()
{
    
	
//	FlT hBath(2.);
	FlT rBath(1.);
	FlT rDisk(.9);
	FlT hDisk(0.1);
	FlT dx(.02);

	FlT rAxis(0.05);
	FlT hAxis(.5);

	FlT wPillar(.2);
	FlT dPillar(.1);

	FlT aCrystal(.5);
	FlT hCrystalBase(.5);
	FlT hCrystalPyramid(.5);

	asl::AVec<int> size(asl::makeAVec(105.,105.,100.));

	asl::AVec<>center(.5*dx*AVec<>(size));

	vector<asl::AVec<>> pillar1{asl::makeAVec(wPillar*.5, dPillar*.5,0.),
	                            asl::makeAVec(-wPillar*.5, dPillar*.5,0.),
	                            asl::makeAVec(-wPillar*.5, -dPillar*.5,0.),
	                            asl::makeAVec(wPillar*.5, -dPillar*.5,0.)};

	vector<asl::AVec<>> pillar2{asl::makeAVec(dPillar*.5, wPillar*.5,0.),
	                            asl::makeAVec(-dPillar*.5, wPillar*.5,0.),
	                            asl::makeAVec(-dPillar*.5, -wPillar*.5,0.),
	                            asl::makeAVec(dPillar*.5, -wPillar*.5,0.)};
	
	vector<asl::AVec<>> pillarC{asl::makeAVec(center[0]+rDisk-dPillar*.5, center[1], 0.),
                                asl::makeAVec(center[0]-rDisk+dPillar*.5, center[1], 0.),
	                            asl::makeAVec(center[0], center[1]+rDisk-dPillar*.5,0.),
	                            asl::makeAVec(center[0], center[1]-rDisk+dPillar*.5,0.)};
	vector<vector<asl::AVec<>>> pillarsPoints(4);
	for(unsigned int i(0); i<4; ++i)
		pillarsPoints[i].resize(4);
	
	for(unsigned int i(0); i<4; ++i)
	{
		pillarsPoints[0][i] = pillar2[i] + pillarC[0];
		pillarsPoints[1][i] = pillar2[i] + pillarC[1];
		pillarsPoints[2][i] = pillar1[i] + pillarC[2];
		pillarsPoints[3][i] = pillar1[i] + pillarC[3];
	}
		

	asl::Block block(size, dx);
	auto mBath(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	auto mPlatform(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	auto mCrystal(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	auto bath(-generateDFCylinderInf(rBath, asl::makeAVec(0.,0.,1.),  
	                                 dx*asl::AVec<FlT>(size)*.5));

	auto diskBottom(generateDFCylinder(rDisk, 
	                                   asl::makeAVec(0., 0., hDisk),  
	                                   asl::makeAVec(center[0], center[1], .5*hDisk)));
	auto diskTop(generateDFCylinder(rDisk, 
	                                asl::makeAVec(0., 0., hDisk),  
	                                asl::makeAVec(center[0], center[1], -.5*hDisk - hAxis + dx*size[2])));
	auto axis(generateDFCylinder(rAxis, 
	                                asl::makeAVec(0., 0., hAxis+hDisk*.5),  
	                                asl::makeAVec(center[0], center[1], - .5*hAxis - hDisk*.25 + dx*size[2])));
	auto dfPillar1(generateDFConvexPolygonPrism(pillarsPoints[0]));
	auto dfPillar2(generateDFConvexPolygonPrism(pillarsPoints[1]));
	auto dfPillar3(generateDFConvexPolygonPrism(pillarsPoints[2]));
	auto dfPillar4(generateDFConvexPolygonPrism(pillarsPoints[3]));
	auto dfPillars((dfPillar1 | dfPillar2 | dfPillar3 | dfPillar4) & 
	               generateDFPlane(makeAVec(0.,0.,-1.), makeAVec(0.,0.,.5*hDisk)) &
	               generateDFPlane(makeAVec(0.,0.,1.), makeAVec(0.,0.,-.5*hDisk - hAxis + dx*size[2])));


	auto crystalB(asl::generateDFConvexPolygonPrism({center+makeAVec( aCrystal,  aCrystal,0.),
						  				             center+makeAVec(-aCrystal,  aCrystal,0.),
										             center+makeAVec(-aCrystal, -aCrystal,0.),
										             center+makeAVec( aCrystal, -aCrystal,0.)}) &
	              generateDFPlane(makeAVec(0.,0.,-1.), makeAVec(0.,0., hDisk)) &
	              generateDFPlane(makeAVec(0.,0., 1.), makeAVec(0.,0., hDisk + hCrystalBase)));
	auto cCrPyrBase(makeAVec(center[0],center[1],hDisk+hCrystalBase-.01));
	auto crystalT(asl::generateDFConvexPolygonPyramid({cCrPyrBase+makeAVec( aCrystal,  aCrystal,0.),
						  				               cCrPyrBase+makeAVec(-aCrystal,  aCrystal,0.),
										               cCrPyrBase+makeAVec(-aCrystal, -aCrystal,0.),
										               cCrPyrBase+makeAVec( aCrystal, -aCrystal,0.)},
	                                                   cCrPyrBase+makeAVec(0.,0.,hCrystalPyramid)));
	
	asl::initData(mBath, normalize(bath, dx));
	asl::initData(mPlatform, normalize(diskBottom | diskTop | axis | dfPillars, dx));
	asl::initData(mCrystal, normalize(crystalB | crystalT, dx));

//	asl::writeVTKXML("distFAdvanced3D.vti", *data, "data");
	asl::WriterVTKXML writer("distFAdvanced3D");
	writer.addScalars("Bath", *mBath);
	writer.addScalars("Platform", *mPlatform);
	writer.addScalars("Crystal", *mCrystal);
	writer.write();

	return true;	
}


bool testDistFOptimizer()
{
	FlT r(10.);
	FlT dx(1.);
	asl::AVec<int> size(asl::makeAVec(50.,50.,50.));

	asl::Block block(size, dx);
	auto data(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	auto center(asl::AVec<FlT>(size)*.5);
	auto df1(generateDFSphere(r, center));
	auto df2(generateDFSphere(r, .6*center));

	asl::initData(data, normalize((df1 | df2),dx));
	optimizeMap(data, &asl::d3q15());

	asl::writeVTKXML("distDistFOptimizer.vti", *data, "data");

	return true;	
}


int main()
{
	testDistFOperations2D();
	testDistFOperations3D();
	testDistFOrderedCylinders();
	testDistFUnorderedCylinders();
	testDistFNormalization2D();
	testDistFNormalization3D();
	testDistFOperations3DPrism();	
	testDistFOperations3DBlock();
	testDistFAdvanced3D();
	testDistFOptimizer();
	
	return 0;
}