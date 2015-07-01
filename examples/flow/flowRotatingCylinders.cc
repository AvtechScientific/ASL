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
	\example flowRotatingCylinders.cc
 */


#include <aslDataInc.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <math/aslPositionFunction.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLBGK.h>
#include <num/aslLBGKBC.h>
#include <utilities/aslTimer.h>
#include <utilities/aslParametersManager.h>


typedef float FlT;
//typedef double FlT;
typedef asl::UValue<double> Param;

using asl::AVec;
using asl::makeAVec;


int main(int argc, char* argv[])
{
	// Optionally add appParamsManager to be able to manipulate at least
	// hardware parameters(platform/device) through command line/parameters file
	asl::ApplicationParametersManager appParamsManager("flowRotatingCylinders",
	                                                   "1.0");
	appParamsManager.load(argc, argv);

	Param dx(1.);
	Param dt(1.);
	Param nu(.01);
	Param w(2e-3);

	Param nuNum(nu.v()*dt.v()/dx.v()/dx.v());
	AVec<int> size(asl::makeAVec(100,100,150));

	AVec<> gSize(dx.v()*AVec<>(size));

	
	std::cout << "Data initialization... ";

	asl::Block block(size,dx.v());

	auto exCylMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(exCylMap, 
	              asl::normalize(-generateDFCylinderInf(.48*gSize[0],
	                                                 makeAVec(0.,0.,1.),
	                                                 .5*gSize), 
	                             dx.v()));
	auto inCylMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(inCylMap, 
	              asl::normalize(generateDFCylinderInf(.24*gSize[0],
	                                                makeAVec(0.,0.,1.),
	                                                .5*gSize),
	                             dx.v()));

	auto computationalDomainMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(computationalDomainMap, 
	              asl::normalize(-generateDFCylinderInf(.48*gSize[0],
	                                                 makeAVec(0.,0.,1.),
	                                                 .5*gSize) | 
	                             generateDFCylinderInf(.24*gSize[0],
	                                                 makeAVec(0.,0.,1.),
	                                                 .5*gSize) |
	                             asl::generateDFInBlock(block, 0), 
	                             dx.v()));

		
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... ";

	asl::SPLBGK lbgk(new asl::LBGK(block, 
				               acl::generateVEConstant(FlT(nuNum.v())),  
	        			       &asl::d3q15()));
	
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(.0,.0,.0));

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcV;

	auto vfEx(asl::generatePFRotationField(makeAVec(0.,0., w.v()), .5*gSize));
	auto vfIn(asl::generatePFRotationField(makeAVec(0.,0.,-2.*w.v()), .5*gSize));	
	
	bc.push_back(generateBCVelocity(lbgk, vfEx, exCylMap,computationalDomainMap));
	bcV.push_back(generateBCNoSlipVel(lbgk, exCylMap));
	bc.push_back(generateBCVelocity(lbgk, vfIn, inCylMap,computationalDomainMap));
	bcV.push_back(generateBCNoSlipVel(lbgk, inCylMap));
	bc.push_back(generateBCNoSlip(lbgk,{asl::Z0, asl::ZE}));

	initAll(bc);
	initAll(bcV);

	std::cout << "Finished" << endl;
	std::cout << "Computing...";
	asl::Timer timer;

	asl::WriterVTKXML writer("flowRotCylRes");
//	writer.addScalars("mapEx", *exCylMap);
//	writer.addScalars("mapIn", *inCylMap);
	writer.addScalars("map", *computationalDomainMap);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *lbgk->getVelocity());

	executeAll(bc);

	executeAll(bcV);
	writer.write();

	timer.start();
	for (unsigned int i(0); i < 10001  ; ++i)
	{
		lbgk->execute();
		executeAll(bc);
		if (!(i%1000))
		{
			cout <<  i  << endl;
			executeAll(bcV);
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
