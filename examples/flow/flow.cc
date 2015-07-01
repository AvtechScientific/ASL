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
	\example flow.cc
 */

#include <utilities/aslParametersManager.h>
#include <aslDataInc.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLBGK.h>
#include <num/aslLBGKBC.h>
#include <utilities/aslTimer.h>
#include <acl/aclUtilities.h>



typedef float FlT;
//typedef double FlT;
typedef asl::UValue<double> Param;

using asl::AVec;
using asl::makeAVec;


int main(int argc, char* argv[])
{
	// Optionally add appParamsManager to be able to manipulate at least
	// hardware parameters(platform/device) through command line/parameters file
	asl::ApplicationParametersManager appParamsManager("flow",
	                                                   "1.0");
	appParamsManager.load(argc, argv);

	Param dx(1.);
	Param dt(1.);
	Param nu(.00625);

	Param nuNum(nu.v()*dt.v()/dx.v()/dx.v());
	AVec<int> size(asl::makeAVec(300, 50, 50));

	auto gSize(dx.v()*AVec<>(size));

	
	std::cout << "Data initialization... ";

	asl::Block block(size,dx.v());

	auto ball(generateDFCylinderInf(.125*gSize[1],
	                             makeAVec(0.,1.,0.),
	                             .45*makeAVec(gSize[1],gSize[1],gSize[2])));
	// Formula
	auto ballMap(asl::generateDataContainer_SP(block, ball, 1u, acl::typeToTypeID<FlT>()));
	// Data
	auto ballMapMem(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(ballMapMem, ball);

	
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... ";

	asl::SPLBGK lbgk(new asl::LBGK(block,
	                               acl::generateVEConstant(FlT(nuNum.v())),
	                               &asl::d3q15()));
	lbgk->init();

	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	// Initialization of lbgk
	lbgkUtil->initF(acl::generateVEConstant(.0,.0,.0));
//	lbgkUtil->initF(acl::generateVEConstant(.1), acl::generateVEConstant(.0,.0,.0));

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcVis;
	
	bc.push_back(asl::generateBCConstantPressure(lbgk, 1.2, {asl::X0}));
	bc.push_back(asl::generateBCConstantPressure(lbgk, 0.8, {asl::XE}));
	bc.push_back(generateBCNoSlip(lbgk, {asl::Y0, asl::YE, asl::Z0, asl::ZE}));
	bc.push_back(generateBCNoSlip(lbgk, ballMap));
	bcVis.push_back(generateBCNoSlipVel(lbgk, ballMap));

	initAll(bc);
	initAll(bcVis);

	std::cout << "Finished" << endl;
	std::cout << "Computing...";
	asl::Timer timer;

	asl::WriterVTKXML writer("flow");
	writer.addScalars("map", *ballMapMem);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *lbgk->getVelocity());



	executeAll(bc);
	executeAll(bcVis);
	writer.write();

	timer.start();
	for (unsigned int i(1); i < 1001; ++i)
	{
		lbgk->execute();
		executeAll(bc);
		if (!(i%5000))
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
