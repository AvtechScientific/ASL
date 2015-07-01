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
	\example flowKDPGrowth.cc
 */

#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <math/aslPositionFunction.h>
#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLBGK.h>
#include <num/aslLBGKBC.h>
#include <num/aslBasicBC.h>
#include <num/aslCrystalGrowthBC.h>
#include <num/aslFDAdvectionDiffusion.h>
#include <utilities/aslTimer.h>

using asl::AVec;
using asl::makeAVec;

asl::SPDistanceFunction generateBath(asl::Block & bl)
{

//	double hBath(2.);
	double rBath(1.);

	double dx(bl.dx);

	asl::AVec<int> size(bl.getSize());

	asl::AVec<>center(.5*dx*AVec<>(size));

		

	auto bath(-(generateDFCylinderInf(rBath, makeAVec(0.,0.,1.), dx*AVec<>(size)*.5) &
	           generateDFPlane(makeAVec(0.,0.,1.), center*1.99) &
	           generateDFPlane(makeAVec(0.,0.,-1.), center*0.)));

	return normalize(bath, dx);
}

asl::SPDistanceFunction generatePlatform(asl::Block & bl)
{
	double rDisk(.9);
	double hDisk(0.1);

	double rAxis(0.05);
	double hAxis(.5);

	double wPillar(.2);
	double dPillar(.1);

	double dx(bl.dx);
	asl::AVec<int> size(bl.getSize());
	asl::AVec<>center(.5*dx*AVec<>(size));
	
	vector<asl::AVec<>> pillar1{makeAVec(wPillar*.5, dPillar*.5,0.),
	                            makeAVec(-wPillar*.5, dPillar*.5,0.),
	                            makeAVec(-wPillar*.5, -dPillar*.5,0.),
	                            makeAVec(wPillar*.5, -dPillar*.5,0.)};

	vector<asl::AVec<>> pillar2{makeAVec(dPillar*.5, wPillar*.5,0.),
	                            makeAVec(-dPillar*.5, wPillar*.5,0.),
	                            makeAVec(-dPillar*.5, -wPillar*.5,0.),
	                            makeAVec(dPillar*.5, -wPillar*.5,0.)};
	
	vector<asl::AVec<>> pillarC{makeAVec(center[0]+rDisk-dPillar*.5, center[1], 0.),
                                makeAVec(center[0]-rDisk+dPillar*.5, center[1], 0.),
	                            makeAVec(center[0], center[1]+rDisk-dPillar*.5,0.),
	                            makeAVec(center[0], center[1]-rDisk+dPillar*.5,0.)};
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

	
	auto diskBottom(generateDFCylinder(rDisk, 
	                                   makeAVec(0., 0., hDisk),  
	                                   makeAVec(center[0], center[1], .5*hDisk)));
	auto diskTop(generateDFCylinder(rDisk, 
	                                makeAVec(0., 0., hDisk),  
	                                makeAVec(center[0], center[1], -.5*hDisk - hAxis + dx*size[2])));
	auto axis(generateDFCylinder(rAxis, 
	                                makeAVec(0., 0., hAxis+hDisk*.5),  
	                                makeAVec(center[0], center[1], - .5*hAxis - hDisk*.25 + dx*size[2])));
	auto dfPillar1(generateDFConvexPolygonPrism(pillarsPoints[0]));
	auto dfPillar2(generateDFConvexPolygonPrism(pillarsPoints[1]));
	auto dfPillar3(generateDFConvexPolygonPrism(pillarsPoints[2]));
	auto dfPillar4(generateDFConvexPolygonPrism(pillarsPoints[3]));
	auto dfPillars((dfPillar1 | dfPillar2 | dfPillar3 | dfPillar4) & 
	               generateDFPlane(makeAVec(0.,0.,-1.), makeAVec(0.,0.,.5*hDisk)) &
	               generateDFPlane(makeAVec(0.,0.,1.), makeAVec(0.,0.,-.5*hDisk - hAxis + dx*size[2])));

	return normalize(diskBottom | diskTop | axis | dfPillars, dx);
}

asl::SPDistanceFunction generateCrystal(asl::Block & bl)
{

	double aCrystal(.5);
	double hCrystalBase(.5);
	double hCrystalPyramid(.5);

	double hDisk(0.1);

	double dx(bl.dx);
	asl::AVec<int> size(bl.getSize());
	asl::AVec<>center(.5*dx*AVec<>(size));
	
	auto crystalB(asl::generateDFConvexPolygonPrism({center+makeAVec( aCrystal,  aCrystal,0.),
						  				             center+makeAVec(-aCrystal,  aCrystal,0.),
										             center+makeAVec(-aCrystal, -aCrystal,0.),
										             center+makeAVec( aCrystal, -aCrystal,0.)}) &
	              generateDFPlane(makeAVec(0.,0.,-1.), makeAVec(0.,0., hDisk-.001)) &
	              generateDFPlane(makeAVec(0.,0., 1.), makeAVec(0.,0., hDisk + hCrystalBase)));
	auto cCrPyrBase(makeAVec(center[0],center[1],hDisk+hCrystalBase-.01));
	auto crystalT(asl::generateDFConvexPolygonPyramid({cCrPyrBase+makeAVec( aCrystal,  aCrystal,0.),
						  				               cCrPyrBase+makeAVec(-aCrystal,  aCrystal,0.),
										               cCrPyrBase+makeAVec(-aCrystal, -aCrystal,0.),
										               cCrPyrBase+makeAVec( aCrystal, -aCrystal,0.)},
	                                                   cCrPyrBase+makeAVec(0.,0.,hCrystalPyramid)));
	return normalize(crystalB | crystalT, dx);
//	return crystalB | crystalT;
}

double getWRotation(double t)
{
	double tPeriod(128);
	double wMax(6.*3.14*2./60.);
	double tPlato(tPeriod * .25);
	double tAcceleration(tPeriod * .1);
	double tStop(tPeriod * .05);

	double intPart;
	double tRel(modf(t/tPeriod, &intPart));
	double x(0);
	if(tRel<=tAcceleration)
		x = tRel / tAcceleration;
	if(tRel>tAcceleration && tRel<=tAcceleration+tPlato)
		x = 1.;
	if(tRel>tAcceleration+tPlato && tRel<=2.*tAcceleration+tPlato)
		x = (2.*tAcceleration + tPlato - tRel) / tAcceleration;
	if(tRel>2.*tAcceleration+tPlato && tRel<=2.*tAcceleration+tPlato+tStop)
		x = 0;
	if(tRel>2.*tAcceleration+tPlato+tStop && tRel<=3.*tAcceleration+tPlato+tStop)
		x = -(tRel-2.*tAcceleration-tPlato-tStop) / tAcceleration;
	if(tRel>3.*tAcceleration+tPlato+tStop && tRel<=3.*tAcceleration+2.*tPlato+tStop)
		x = -1.;
	if(tRel>3.*tAcceleration+2.*tPlato+tStop && tRel<=4.*tAcceleration+2.*tPlato+tStop)
		x = -(4.*tAcceleration+2.*tPlato+tStop-tRel)/tAcceleration;
	if(tRel>4.*tAcceleration+2.*tPlato+tStop)
		x = 0;
	return wMax*x;

//	flux = -9.32e-5*(1.170 - c); c_0=0.326 ceq=0.267
}


typedef float FlT;
//typedef double FlT;
typedef asl::UValue<double> Param;

using asl::AVec;
using asl::makeAVec;


int main(int argc, char* argv[])
{
	// Optionally add appParamsManager to be able to manipulate at least
	// hardware parameters(platform/device) through command line/parameters file
	asl::ApplicationParametersManager appParamsManager("flowKDPGrowth",
	                                                   "1.0");
	appParamsManager.load(argc, argv);

	Param dx(.02);
	Param dt(0.8e-2);
	Param nu(1e-2);
	Param difC(1e-2/300.);
//	Param w(48.*3.14*2./60.);
	// Angular velocity
	Param w(6.*3.14*2./60.);
	
	Param nuNum(nu.v()*dt.v()/dx.v()/dx.v());
	Param difCNum(difC.v()*dt.v()/dx.v()/dx.v());

	// Angular velocity in one iteration
	Param wNum(w.v()*dt.v());

	Param c0(0.326);
	
	AVec<int> size(asl::makeAVec(105.,105.,100.));

	AVec<> gSize(dx.v()*AVec<>(size));

	std::cout << "Data initialization...";

    auto templ(&asl::d3q19());	
	asl::Block block(size,dx.v());

	auto bathMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(bathMap, generateBath(block));
	auto platformCrysMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(platformCrysMap, generatePlatform(block) | generateCrystal(block));
	auto bathPlatformMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(bathPlatformMap, generateBath(block) | generatePlatform(block));
	auto bathPlatformCrystalMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(bathPlatformCrystalMap, generateBath(block) | generatePlatform(block) | generateCrystal(block));
	auto crystalMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(crystalMap, generateCrystal(block));

	auto cField(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(cField, c0.v());
	
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization...";

	asl::SPLBGK lbgk(new asl::LBGK(block, 
				               acl::generateVEConstant(FlT(nuNum.v())),  
	        			       templ));
	// Set angular velocity in lbgk
	lbgk->setOmega(acl::generateVEConstant(makeAVec(0.,0.,wNum.v())));
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(.0,.0,.0));

	auto nmDif(asl::generateFDAdvectionDiffusion(cField, 
	                                             difCNum.v(), 
	                                             lbgk->getVelocity(), 
	                                             templ,
	                                             true));

	nmDif->init();
	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcV;
	std::vector<asl::SPNumMethod> bcDif;

	// Position Function Angular Velocity Field
	auto vfBath(asl::generatePFRotationField(makeAVec(0.,0., wNum.v()/dx.v()), .5*gSize));
	// Boundary condition
	bc.push_back(generateBCVelocity(lbgk, vfBath, bathMap));
	// Boundary condition for visualization
	bcV.push_back(generateBCVelocityVel(lbgk, vfBath, bathMap));
	bc.push_back(generateBCNoSlip(lbgk,  platformCrysMap));
	bcV.push_back(generateBCNoSlipVel(lbgk, platformCrysMap));
	bcDif.push_back(generateBCConstantGradient(cField, 
	                                        0., 
	                                        bathPlatformMap, 
	                                        bathPlatformCrystalMap, 
	                                        templ));
	bcDif.push_back(generateBCLinearGrowth2(cField, 1.17, 
	                                        -9.32e-6/difC.v()*dx.v(), 
	                                        crystalMap, 
	                                        bathPlatformCrystalMap, 
	                                        templ));
//	bcDif.push_back(generateBCConstantGradient2(cField, .1, crystalMap, bathPlatformCrystalMap, templ));
	initAll(bc);
	initAll(bcV);
	initAll(bcDif);

	std::cout << "Finished" << endl;
	std::cout << "Computing...";
	asl::Timer timer;
	asl::Timer timerBC;

	asl::WriterVTKXML writer("flowKDPGrowthRes");
	writer.addScalars("mapBath", *bathMap);
	writer.addScalars("mapPlatformCrys", *platformCrysMap);
	writer.addScalars("mapBathPlatformCrystal", *bathPlatformCrystalMap);
	writer.addScalars("mapCrys", *crystalMap);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addScalars("c", *cField);
	writer.addVector("v", *lbgk->getVelocity());

	executeAll(bcV);
	executeAll(bc);
	executeAll(bcDif);

	writer.write();

	timer.start();
	timerBC.reset();
	for (unsigned int i(0); i <= 8001  ; ++i)
	{
		lbgk->execute();
		timerBC.resume();
		executeAll(bcV);
		executeAll(bc);
		timerBC.stop();
		nmDif->execute();
		timerBC.resume();
		executeAll(bcDif);
		timerBC.stop();
		
		if (!(i%2000))
		{
			cout <<  i  << endl;
			writer.write();
		}
	}
	timer.stop();
	
	std::cout << "Finished" << endl;	

	cout << "time=" << timer.getTime() << "; clockTime="
		 <<  timer.getClockTime() <<  "; load=" 
		 <<  timer.getProcessorLoad() * 100 << "%; timeBC = " 
		 <<  timerBC.getTime() << endl;

	std::cout << "Output...";
	std::cout << "Finished" << endl;	
	std::cout << "Ok" << endl;

	return 0;
}
