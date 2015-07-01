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
	\example flow2.cc
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

asl::SPDistanceFunction generateMirror(double x, double y)
{
	double hCyl (10.);
	double rCyl (1.5);
	double a(10.);
	double b(6.);
	vector<AVec<>> rect({.5*makeAVec( .866*a,-b,-a*.5),.5*makeAVec(-.866*a,-b, a*.5),
	                     .5*makeAVec(-.866*a, b, a*.5),.5*makeAVec( .866*a, b,-a*.5)});
	vector<AVec<>> r1(rect.size());
	vector<AVec<>> r2(rect.size());
	auto r1C(makeAVec(x,y+rCyl*.5+b*.5,hCyl));
	auto r2C(makeAVec(x,y-rCyl*.5-b*.5,hCyl));
	for(unsigned int i(0); i<rect.size();++i)
	{
		r1[i]=r1C + rect[i];
		r2[i]=r2C + rect[i];
	}
	
	auto res(asl::generateDFConvexPolygonPrism(r1) | generateDFConvexPolygonPrism(r2));
	res = (res & 
		   generateDFPlane(makeAVec(.5,0.,.866), r1C + makeAVec(rCyl,0.,0.)) &
		   generateDFPlane(makeAVec(-.5,0.,-.866), r1C - makeAVec(rCyl,0.,0.))) |
		   generateDFCylinder(rCyl, makeAVec(0.,0.,hCyl), makeAVec(x,y,hCyl*.5));
	return res;
}

asl::SPDistanceFunction generateMirrors()
{
	vector<double> xValues;
	vector<double> yValues;

	for(unsigned int i(0); i < 8; ++i)
		for(unsigned int j(0); j < 10; ++j)
		{
			xValues.push_back(25.+ 20.*i);
			yValues.push_back(25.+ 15.*j);
		}

	auto res(generateMirror(xValues[0], yValues[0]));
	for(unsigned int i(1); i < xValues.size(); ++i)
		res = res | generateMirror(xValues[i], yValues[i]);

	return res;
}


int main(int argc, char* argv[])
{
	// Optionally add appParamsManager to be able to manipulate at least
	// hardware parameters(platform/device) through command line/parameters file
	asl::ApplicationParametersManager appParamsManager("flow2",
	                                                   "1.0");
	appParamsManager.load(argc, argv);

	Param dx(1.);
	Param dt(1.);
	Param nu(.0125);

	Param nuNum(nu.v()*dt.v()/dx.v()/dx.v());
	AVec<int> size(asl::makeAVec(200,200,25));

	AVec<> gSize(dx.v()*AVec<>(size));

	
	std::cout << "Data initialization... ";

	asl::Block block(size,dx.v());

	auto mirrorsMapMem(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(mirrorsMapMem, generateMirrors());

	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... ";

	asl::SPLBGK lbgk(new asl::LBGK(block, 
				               acl::generateVEConstant(FlT(nuNum.v())),  
	        			       &asl::d3q15()));
	
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(.0,.0,.0));
//	lbgkUtil->initF(acl::generateVEConstant(.1), acl::generateVEConstant(.0,.0,.0));

	auto bcNoSlip(generateBCNoSlip(lbgk,{asl::Y0, asl::YE, asl::Z0}));
	auto bcNoSlipM(generateBCNoSlip(lbgk, mirrorsMapMem));
	auto bcNoSlipV(generateBCNoSlipVel(lbgk, mirrorsMapMem));
	auto bcP(asl::generateBCConstantPressure(lbgk, 1., {asl::X0, asl::XE}));
	auto bcTop(asl::generateBCConstantVelocity(lbgk, asl::makeAVec(0.1,0.,0.), {asl::ZE}));

	bcNoSlip->init();
	bcNoSlipM->init();
	bcNoSlipV->init();
	bcP->init();
	bcTop->init();

	std::cout << "Finished" << endl;
	std::cout << "Computing...";
	asl::Timer timer;

	asl::WriterVTKXML writer("flow2");
	writer.addScalars("map", *mirrorsMapMem);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *lbgk->getVelocity());



	bcP->execute();
	bcTop->execute();
	bcNoSlip->execute();
	bcNoSlipM->execute();

	bcNoSlipV->execute();
	writer.write();

	timer.start();
	for (unsigned int i(0); i < 1000  ; ++i)
	{
		lbgk->execute();
		bcP->execute();
		bcTop->execute();
		bcNoSlip->execute();
		bcNoSlipM->execute();
		if (!(i%100))
		{
			cout <<  i  << endl;
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
