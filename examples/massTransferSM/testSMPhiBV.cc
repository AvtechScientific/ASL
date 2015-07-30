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
	\example testSMPhiBV.cc
 */

#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslFDStefanMaxwell.h>
#include <num/aslFDAdvectionDiffusion.h>
#include <num/aslFDElChemKinetics.h>
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
	asl::ApplicationParametersManager appParamsManager("testSMPhiBV",
	                                                   "1.0");
	appParamsManager.load(argc, argv);

	Param dx(1.);
	Param dt(1.);
	Param diffCoef(.15);
	
	Param diffCoefNum(diffCoef.v()*dt.v()/dx.v()/dx.v());
	
	AVec<int> size(asl::makeAVec(10,20,20));

	auto gSize(dx.v()*AVec<>(size));

	std::cout << "Data initialization... ";

	asl::Block block(size,dx.v());
	
	auto c1Field(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(c1Field, 0.5);
	auto c2Field(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(c2Field, 0.5);	
	auto c1aField(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(c1aField, 0.5);	

	auto phiField(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(phiField, 0.5);	

	auto phiSField(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(phiSField, 0.);	
	
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... ";

	auto templ(&asl::d3q7());
	auto nm(generateFDStefanMaxwell(c1Field, c2Field,  diffCoefNum.v(), templ));
	nm->setDustDiffusionCoefficient(0,acl::generateVEConstant(diffCoefNum.v()*.5));
	nm->setDustDiffusionCoefficient(1,acl::generateVEConstant(diffCoefNum.v()));
	nm->setElectricField(asl::generateDataContainer_SP(block,
	                                                   phiField->getEContainer()*1e5/8.31/300,
	                                                   1u));
	nm->setCharge(0,acl::generateVEConstant(-2.));
	nm->setCharge(1,acl::generateVEConstant(0.));
	nm->init();

	auto nmA(generateFDAdvectionDiffusion(c1aField, diffCoefNum.v(), templ));
	nmA->init();
	
	auto nmPhi(make_shared<asl::FDStefanMaxwellElectricField>(nm, phiField));
	nmPhi->setPhiS(asl::generateDataContainer_SP(block,
	                                             -phiSField->getEContainer()*.1,
	                                             1u));
	nmPhi->init();

	auto nmBV(generateFDBVKinetics(c1aField, -1., c1Field, 1., phiSField, 
	                               asl::generateDataContainer_SP(block,
	                                                             phiField->getEContainer()+.2,
	                                                             1u),
	                               1., .5, 1.));
	nmBV->addAI(c2Field, -1);
	nmBV->init();
	
	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcPhi;	

	bc.push_back(asl::generateBCConstantGradient(c1Field, 0, templ, {asl::XE,asl::Y0,asl::YE,asl::Z0,asl::ZE}));
	bc.push_back(asl::generateBCConstantValue(c1Field, 0., {asl::X0}));
	bc.push_back(asl::generateBCConstantGradient(c2Field, 0, templ, {asl::XE,asl::Y0,asl::YE,asl::Z0,asl::ZE}));
	bc.push_back(asl::generateBCConstantValue(c2Field, 1, {asl::X0}));
	bc.push_back(asl::generateBCConstantGradient(c1aField, 0, templ, {asl::X0,asl::Y0,asl::YE,asl::Z0,asl::ZE}));
	bc.push_back(asl::generateBCConstantValue(c1aField, 1, {asl::XE}));
	initAll(bc);
	bcPhi.push_back(asl::generateBCConstantGradient(phiField, 0, templ, {asl::XE,asl::Y0,asl::YE,asl::Z0,asl::ZE}));
	bcPhi.push_back(asl::generateBCConstantValue(phiField, 0, {asl::X0}));
	initAll(bcPhi);

	std::cout << "Finished" << endl;
	std::cout << "Computing..." << flush;
	asl::Timer timer;

	asl::WriterVTKXML writer(appParamsManager.getDir() + "testSMPhiBV");
	writer.addScalars("c1", *c1Field);
	writer.addScalars("c2", *c2Field);
	writer.addScalars("c1a", *c1aField);	
	writer.addScalars("phi", *phiField);
	writer.addScalars("phiS", *phiSField);

	executeAll(bc);
	executeAll(bcPhi);
	writer.write();

	timer.start();
	for (unsigned int i(1); i < 401; ++i)
	{
		for (unsigned int j(0); j<50; ++j)
		{
			nmBV->executeJ();
			nmPhi->execute();
			executeAll(bcPhi);		
		}
		nm->execute();
		nmA->execute();
		nmBV->execute();
		executeAll(bc);
		if (!(i%40))
		{
			cout << i << endl;
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
