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
	\example flow3.cc
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

asl::SPDistanceFunction generateOrderedCylinders(asl::Block & block)
{
	double r(2.);
	double spacing(4.);

	asl::SPDistanceFunction cylinder;
	asl::SPDistanceFunction resultGeometry;
	asl::AVec<double> orientation(asl::makeAVec(0., 0., 1.));
	for (int i = 0; i < block.getSize()[0] / (2 * r + spacing); ++i)
	{
		for (int j = 0; j < block.getSize()[1] / (2 * r + spacing); ++j)
		{
			cylinder = generateDFCylinderInf(r, orientation, asl::makeAVec(i * (2. * r + spacing) + r + spacing / 2., j * (2. * r + spacing) + r + spacing / 2., 0.));
			resultGeometry = resultGeometry | cylinder;
		}
	}

	return resultGeometry;
}


int main(int argc, char* argv[])
{
	// Optionally add appParamsManager to be able to manipulate at least
	// hardware parameters(platform/device) through command line/parameters file
	asl::ApplicationParametersManager appParamsManager("flow3",
	                                                   "1.0");
	appParamsManager.load(argc, argv);

	Param dx(1.);
	Param dt(1.);
	Param nu(.00625);

	Param nuNum(nu.v()*dt.v()/dx.v()/dx.v());
	AVec<int> size(asl::makeAVec(50, 50, 50));

	auto gSize(dx.v()*AVec<>(size));

	
	std::cout << "Flow: Data initialization...";

	asl::Block block(size,dx.v());

	auto cylindersMapMem(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(cylindersMapMem, generateOrderedCylinders(block));

	
	std::cout << "Finished" << endl;
	
	std::cout << "Flow: Numerics initialization...";

	asl::SPLBGK lbgk(new asl::LBGK(block, 
				               acl::generateVEConstant(FlT(nuNum.v())),  
	        			       &asl::d3q15()));
	
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(.0,.0,.0));

	auto bcNoSlip(generateBCNoSlip(lbgk,{asl::Y0, asl::YE, asl::Z0, asl::ZE}));
	auto bcNoSlipM(generateBCNoSlip(lbgk, cylindersMapMem));
	auto bcNoSlipV(generateBCNoSlipVel(lbgk, cylindersMapMem));
	asl::BCConstantPressure bcIn(lbgk, acl::generateVEConstant(1.2));
	asl::BCConstantPressure bcOut(lbgk, acl::generateVEConstant(0.8));

	addSliceX0(bcIn);
	addSliceXE(bcOut);

	bcNoSlip->init();
	bcNoSlipM->init();
	bcNoSlipV->init();
	bcIn.init();
	bcOut.init();

	std::cout << "Finished" << endl;
	std::cout << "Computing...";
	asl::Timer timer;

	asl::WriterVTKXML writer("flow3Res");
	writer.addScalars("map", *cylindersMapMem);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *lbgk->getVelocity());



	bcIn.execute();
	bcOut.execute();
	bcNoSlip->execute();
	bcNoSlipM->execute();

	bcNoSlipV->execute();
	writer.write();

	timer.start();
	for (unsigned int i(0); i < 1000; ++i)
	{
		lbgk->execute();
		bcIn.execute();
		bcOut.execute();
		bcNoSlip->execute();
		bcNoSlipM->execute();
		if (!(i%100))
		{
			cout << i << endl;
			bcNoSlipV->execute();
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
