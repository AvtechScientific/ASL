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
	\example locomotive_stability.cc
	Required input file: [locomotive.stl](http://asl.org.il/input_data/locomotive.stl)
 */

#include <math/aslVectors.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <math/aslPositionFunction.h>
#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLBGK.h>
#include <num/aslLBGKBC.h>
#include <utilities/aslTimer.h>
#include <acl/aclUtilities.h>
#include <readers/aslVTKFormatReaders.h>
#include <utilities/aslParametersManager.h>


typedef float FlT;
//typedef double FlT;
typedef asl::UValue<double> Param;

using asl::AVec;
using asl::makeAVec;

asl::SPDistanceFunction generateNozzels(asl::Block & bl)
{

	double y(10.33);
	double z(3.08);
	double h(1.);
	double rNozzel(.3);
	double n(4.);
	double x0(-20.);
	double xE(-6.);	

	double dx(bl.dx);

	asl::SPDistanceFunction res(generateDFCylinder(rNozzel, 
	                                               makeAVec(0., 0., h), 
	                                               makeAVec(x0, y, z + h * .5)));
	for(unsigned int i(1); i < n; ++i)
		res = res | generateDFCylinder(rNozzel, 
	                                   makeAVec(0., 0., h), 
	                                   makeAVec(x0 + i * (xE-x0)/n, y, z+h*.5)); 
	return normalize(res, dx);
}


int main(int argc, char* argv[])
{
	// Optionally add appParamsManager to be able to manipulate at least
	// hardware parameters(platform/device) through command line/parameters file
	asl::ApplicationParametersManager appParamsManager("locomotive_stability",
	                                                   "1.0");
	asl::Parameter<string> input("input", "path to the geometry input file");
	appParamsManager.load(argc, argv);
	
	Param dx(0.25);
	Param dt(1.);
	Param nu(.01);
	
	Param nuNum(nu.v()*dt.v()/dx.v()/dx.v());
	
	std::cout << "Data initialization... " << flush;

	auto locomotive(asl::readSurface(input.v(), dx.v(), .5, 1., 0., 1., 2., 4.));
	
	asl::Block block(locomotive->getInternalBlock());

	auto forceField(asl::generateDataContainerACL_SP<FlT>(block, 3, 1u));
	asl::initData(forceField, makeAVec(0., 0., 0.));
	
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... " << flush;

	asl::SPLBGK lbgk(new asl::LBGK(block, 
				               acl::generateVEConstant(FlT(nu.v())),  
	        			       &asl::d3q15()));
	
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(.1,.0,.05));

	auto vfTunnel(asl::generatePFConstant(makeAVec(0.1,0.,0.)));

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcV;

	auto nozzelsMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(nozzelsMap, generateNozzels(block));

	
	bc.push_back(generateBCNoSlip(lbgk,  locomotive));
	bcV.push_back(generateBCNoSlipVel(lbgk, locomotive));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1., 
	                                                makeAVec(0.1,0.,0.05), 
	                                                {asl::X0, asl::XE,asl::Y0,asl::YE,asl::Z0,asl::ZE}));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1., makeAVec(0.,0.,-0.1),  nozzelsMap));
	initAll(bc);
	initAll(bcV);

	auto computeForce(generateComputeSurfaceForce(lbgk, forceField, locomotive));
	computeForce->init();
	

	std::cout << "Finished" << endl;
	std::cout << "Computing..." << endl;

	asl::WriterVTKXML writer(appParamsManager.getDir() + "locomotive_stability");
	writer.addScalars("locomotive", *locomotive);
	writer.addScalars("nozzels", *nozzelsMap);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *lbgk->getVelocity());
	writer.addVector("force", *forceField);

	executeAll(bc);
	executeAll(bcV);
	computeForce->execute();

	writer.write();

	asl::Timer timer, timer1, timer2;
	timer.start();
	for(unsigned int i(1); i < 40001; ++i)
	{
		lbgk->execute();
		executeAll(bc);
		if(!(i%1000))
		{
			cout << i << endl;
			executeAll(bcV);
			computeForce->execute();
			writer.write();
		}
	}
	timer.stop();
	
	std::cout << "Finished" << endl;	

	cout << "time=" << timer.getTime() << "; clockTime="
		 <<  timer.getClockTime()	 <<  "; load=" 
		 <<  timer.getProcessorLoad() * 100 << "%" << endl;

	std::cout << "Output...";
	std::cout << "Finished" << endl;	
	std::cout << "Ok" << endl;

	return 0;
}