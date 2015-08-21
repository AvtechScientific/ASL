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
	\example bus_wind.cc
	Required input file: [bus.stl](http://asl.org.il/input_data/bus.stl)
 */

#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <math/aslDistanceFunction.h>
#include <math/aslPositionFunction.h>
#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLBGK.h>
#include <num/aslLBGKBC.h>
#include <utilities/aslTimer.h>
#include <readers/aslVTKFormatReaders.h>



typedef float FlT;
//typedef double FlT;
typedef asl::UValue<double> Param;

using asl::AVec;
using asl::makeAVec;


int main(int argc, char* argv[])
{
	// Optionally add appParamsManager to be able to manipulate at least
	// hardware parameters(platform/device) through command line/parameters file
	asl::ApplicationParametersManager appParamsManager("bus_wind",
	                                                   "1.0");
	asl::Parameter<string> input("input", "path to the bus geometry input file");
	appParamsManager.load(argc, argv);

	Param dx(8);
	Param dt(1.);
	Param nu(.01);
	
	Param nuNum(nu.v()*dt.v()/dx.v()/dx.v());
	
	std::cout << "Data initialization... " << flush;


	auto object(asl::readSurface(input.v(), dx.v(), 1.5,.25,0.,1.,3.,1.));
	
	asl::Block block(object->getInternalBlock());

	auto forceField(asl::generateDataContainerACL_SP<FlT>(block, 3, 1u));
	asl::initData(forceField, makeAVec(0.,0.,0.));
	
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... " << flush;

	asl::SPLBGK lbgk(new asl::LBGK(block, 
				               acl::generateVEConstant(FlT(nu.v())),  
	        			       &asl::d3q15()));
	
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(-.1,.0,-.05));

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcV;

	bc.push_back(generateBCNoSlip(lbgk,  object));
	bcV.push_back(generateBCNoSlipVel(lbgk, object));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1., 
	                                                makeAVec(-0.1,0.,-0.05), 
	                                                {asl::X0, asl::XE,asl::Y0,asl::YE,asl::Z0,asl::ZE}));
	initAll(bc);
	initAll(bcV);

	auto computeForce(generateComputeSurfaceForce(lbgk, forceField, object));
	computeForce->init();
	

	std::cout << "Finished" << endl;
	std::cout << "Computing..." << endl;

	asl::WriterVTKXML writer("bus_wind");
	writer.addScalars("bus", *object);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *lbgk->getVelocity());
	writer.addVector("force", *forceField);

	executeAll(bc);
	executeAll(bcV);
	computeForce->execute();

	writer.write();

	asl::Timer timer, timer1, timer2;
	timer.start();
	timer1.reset();
	timer2.reset();
	for (unsigned int i(1); i < 101; ++i)
	{
		timer1.resume();
		lbgk->execute();
		timer1.stop();
		timer2.resume();
		executeAll(bc);
		timer2.stop();
		if (!(i%1000))
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
		 <<  timer.getClockTime() <<  "; load=" 
		 <<  timer.getProcessorLoad() * 100 << "%" << endl;
	cout << "time1=" << timer1.getTime() << "; time2=" << timer2.getTime() << endl;

	std::cout << "Output...";
	std::cout << "Finished" << endl;	
	std::cout << "Ok" << endl;

	return 0;
}
