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
	\example testSMDiff3C.cc
	ternary system
 */

#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslFDStefanMaxwell.h>
#include <num/aslBasicBC.h>
#include <utilities/aslTimer.h>
#include <acl/aclMath/aclVectorOfElements.h>

typedef float FlT;
//typedef double FlT;
typedef asl::UValue<double> Param;

using asl::AVec;
using asl::makeAVec;


int main(int argc, char* argv[])
{
	// Optionally add appParamsManager to be able to manipulate at least
	// hardware parameters(platform/device) through command line/parameters file
	asl::ApplicationParametersManager appParamsManager("testSMDiff3C",
	                                                   "1.0");
	appParamsManager.load(argc, argv);

	Param dx(1.);
	Param dt(1.);
	Param diffCoef(.15);
	
	Param diffCoefNum(diffCoef.v()*dt.v()/dx.v()/dx.v());
	
	AVec<int> size(asl::makeAVec(10,20,20));

	auto gSize(dx.v()*AVec<>(size));

	
	std::cout << "Flow: Data initialization... ";

	asl::Block block(size,dx.v());
	
	auto c1Field(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(c1Field, 0.33);
	auto c2Field(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(c2Field, 0.33);	
	auto c3Field(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(c3Field, 0.33);	

	
	std::cout << "Finished" << endl;
	
	std::cout << "Flow: Numerics initialization... ";

	auto templ(&asl::d3q7());
	auto nm(generateFDStefanMaxwell(c1Field, c2Field,  diffCoefNum.v(), templ));
	nm->addComponent(c3Field,acl::generateVEConstant(asl::makeAVec(diffCoefNum.v()*.5,
	                                                               diffCoefNum.v())));
	nm->setDustDiffusionCoefficient(0,acl::generateVEConstant(diffCoefNum.v()*.5));
	nm->setDustDiffusionCoefficient(1,acl::generateVEConstant(diffCoefNum.v()*.7));
	nm->setDustDiffusionCoefficient(2,acl::generateVEConstant(diffCoefNum.v()));
	nm->init();

	std::vector<asl::SPNumMethod> bc;

	bc.push_back(asl::generateBCConstantGradient(c1Field, 0, templ, {asl::Y0,asl::YE,asl::Z0,asl::ZE}));
	bc.push_back(asl::generateBCConstantValue(c1Field, 0, {asl::X0}));
	bc.push_back(asl::generateBCConstantValue(c1Field, 1, {asl::XE}));
	bc.push_back(asl::generateBCConstantGradient(c2Field, 0, templ, {asl::Y0,asl::YE,asl::Z0,asl::ZE}));
	bc.push_back(asl::generateBCConstantValue(c2Field, 0, {asl::XE}));
	bc.push_back(asl::generateBCConstantValue(c2Field, 1, {asl::X0}));
	bc.push_back(asl::generateBCConstantGradient(c3Field, 0, templ, {asl::X0,asl::XE,asl::Z0,asl::ZE}));
	bc.push_back(asl::generateBCConstantValue(c3Field, 0, {asl::YE}));
	bc.push_back(asl::generateBCConstantValue(c3Field, 1, {asl::Y0}));
	initAll(bc);

	std::cout << "Finished" << endl;
	std::cout << "Computing..." << flush;
	asl::Timer timer;

	asl::WriterVTKXML writer(appParamsManager.getDir() + "testSMDiff3C");
	writer.addScalars("c1", *c1Field);
	writer.addScalars("c2", *c2Field);
	writer.addScalars("c3", *c3Field);

	executeAll(bc);
	writer.write();

	timer.start();
	for(unsigned int i(1); i < 401; ++i)
	{
		nm->execute();
		executeAll(bc);
		if(!(i%40))
		{
			cout << i << endl;
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
