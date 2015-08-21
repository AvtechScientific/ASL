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
	\example compressor.cc
	Required input file: [axial-compressor.stl](http://asl.org.il/input_data/axial-compressor.stl)
 */

#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <aslDataInc.h>
#include <math/aslPositionFunction.h>
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

asl::SPDistanceFunction generateCase(asl::Block & bl)
{
	double rCase((bl.getBPosition()[1] - bl.position[1]) / 2.);

	asl::AVec<> center(.5*(bl.getBPosition() + bl.position)); 
	center[1] = bl.position[1]  + rCase;

	auto comprCase(-(generateDFCylinderInf(rCase, makeAVec(0., 0., 1.), center)));

	return normalize(comprCase, bl.dx);
}


int main(int argc, char* argv[])
{
	// Optionally add appParamsManager to be able to manipulate at least
	// hardware parameters(platform/device) through command line/parameters file
	asl::ApplicationParametersManager appParamsManager("compressor",
	                                                   "1.0");
	asl::Parameter<string> input("input", "path to the compressor geometry input file");
	appParamsManager.load(argc, argv);

	Param dx(0.5);
	Param dt(.001);
	Param nu(.2);
	// Angular velocity
	Param w(6.*3.14*2./60.);
	
	AVec<int> size(makeAVec(33., 33., 90.)*(1./dx.v()));
	asl::Block bl(size, dx.v(), makeAVec(-16.65, -26.15, -67.7));
	
//	Param nuNum(nu.v()*dt.v()/dx.v()/dx.v());
	Param nuNum(.2);
	// Angular velocity in one iteration
	Param wNum(w.v()*dt.v());

	std::cout << "Data initialization..." << flush;


	auto compressorMap(asl::readSurface(input.v(), bl));
	
	asl::Block block(compressorMap->getInternalBlock());

	auto comprCaseMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(comprCaseMap, generateCase(block));

	
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization..." << flush;

	asl::SPLBGK lbgk(new asl::LBGK(block,
	                               acl::generateVEConstant(FlT(nuNum.v())),
	                               &asl::d3q15()));
	// Set angular velocity in lbgk
	lbgk->setOmega(acl::generateVEConstant(makeAVec(0.,0.,wNum.v())));
	lbgk->init();

	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	// Initialization of lbgk
	lbgkUtil->initF(acl::generateVEConstant(.0,.0,.0));

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcVis;

	// Position Function Angular Velocity Field
	auto vfCase(asl::generatePFRotationField(makeAVec(0.,0., wNum.v()/dx.v()),
	                                         .5*(block.getBPosition() + block.position)));
	// Boundary condition
	bc.push_back(generateBCVelocity(lbgk, vfCase, comprCaseMap));
	// Boundary condition for visualization
	bcVis.push_back(generateBCVelocityVel(lbgk, vfCase, comprCaseMap));
	bc.push_back(asl::generateBCConstantPressure(lbgk, 1., {asl::Z0, asl::ZE}));
	bc.push_back(generateBCNoSlip(lbgk, compressorMap));
	bcVis.push_back(generateBCNoSlipVel(lbgk, compressorMap));

	initAll(bc);
	initAll(bcVis);


	std::cout << "Finished" << endl;
	std::cout << "Computing..." << endl;
	asl::Timer timer;

	asl::WriterVTKXML writer("compressor");
	writer.addScalars("compressor", *compressorMap);
	writer.addScalars("case", *comprCaseMap);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *lbgk->getVelocity());

	executeAll(bc);
	executeAll(bcVis);

	writer.write();

	timer.start();
	for (unsigned int i(1); i < 10001; ++i)
	{
		lbgk->execute();
		executeAll(bc);
		if (!(i%2000))
		{
			cout << i << endl;
			executeAll(bcVis);
			writer.write();
		}
	}
	timer.stop();
	
	std::cout << "Finished" << endl;	

	cout << "time=" << timer.getTime() << "; clockTime="
		 <<  timer.getClockTime() <<  "; load=" 
		 <<  timer.getProcessorLoad() * 100 << "%" << endl;

	std::cout << "Output...";
	std::cout << "Finished" << endl;	
	std::cout << "Ok" << endl;

	return 0;
}