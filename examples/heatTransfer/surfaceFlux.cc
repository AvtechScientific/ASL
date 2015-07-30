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
	\example surfaceFlux.cc
 */

#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslFDAdvectionDiffusion.h>
#include <num/aslBasicBC.h>
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
	asl::ApplicationParametersManager appParamsManager("surfaceFlux",
	                                                   "1.0");
	appParamsManager.load(argc, argv);

	Param dx(1.);
	Param dt(1.);
	Param diffCoef(.15);
	
	Param diffCoefNum(diffCoef.v()*dt.v()/dx.v()/dx.v());
	AVec<int> size(asl::makeAVec(50,50,50));

	auto gSize(dx.v()*AVec<>(size));

	
	std::cout << "Data initialization... ";

	asl::Block block(size,dx.v());

	auto ball(normalize(//generateDFPlane(makeAVec(1.,3.,1.),.5*gSize+makeAVec(.2,.0,.6)) & 
	                    generateDFSphere(10, .5*gSize+makeAVec(.2,.1,.3)),
	                    dx.v()));
	// Formula
//	auto ballMap(asl::generateDataContainer_SP(block, ball, 1u, acl::typeToTypeID<FlT>()));
	// Data
	auto ballMapMem(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(ballMapMem, ball);

	auto ballB(normalize(-generateDFSphere(24, .5*gSize),dx.v() ));
	auto ballBMapMem(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(ballBMapMem, ballB);
	
	auto cField(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(cField, 0.);

	
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... ";

	auto templ(&asl::d3q15());
	auto nm(generateFDAdvectionDiffusion(cField,  diffCoefNum.v(), templ));
	nm->init();

	std::vector<asl::SPNumMethod> bc;
	
//	bc.push_back(generateBCConstantGradient(cField, 0.1, ballMapMem, templ));
	bc.push_back(generateBCConstantGradient2(cField, 0.1, ballMapMem, templ));
	bc.push_back(asl::generateBCConstantValue(cField, 0, ballBMapMem));
	initAll(bc);

	std::cout << "Finished" << endl;
	std::cout << "Computing..." << flush;
	asl::Timer timer;

	asl::WriterVTKXML writer(appParamsManager.getDir() + "surfaceFlux");
	writer.addScalars("map", *ballMapMem);
	writer.addScalars("mapE", *ballBMapMem);
	writer.addScalars("c", *cField);

	executeAll(bc);
	writer.write();

	timer.start();
	for (unsigned int i(1); i < 201; ++i)
	{
		nm->execute();
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
		 <<  timer.getClockTime() << "; load=" 
		 <<  timer.getProcessorLoad() * 100 << "%" << endl;

	std::cout << "Output...";
	std::cout << "Finished" << endl;	
	std::cout << "Ok" << endl;

	return 0;
}
