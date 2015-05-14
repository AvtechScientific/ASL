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

#include <utilities/aslUValue.h>
#include <math/aslVectors.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <math/aslPositionFunction.h>
#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLBGK.h>
#include <num/aslLBGKBC.h>
#include "utilities/aslTimer.h"
#include "acl/aclUtilities.h"
#include <readers/aslVTKFormatReaders.h>



typedef float FlT;
//typedef double FlT;
typedef asl::UValue<double> Param;

using asl::AVec;
using asl::makeAVec;

asl::SPDistanceFunction generateTunel(asl::Block & bl)
{

	double l(bl.getBPosition()[0]-bl.position[0]+bl.dx);
	double rTunel((bl.getBPosition()[2]-bl.position[2])/2.);

	double dx(bl.dx);

	asl::AVec<int> size(bl.getSize());

	asl::AVec<> center(.5*(bl.getBPosition()+bl.position)); 
	center[1]=bl.position[1]+.25*rTunel;
	asl::AVec<> centerG(center); 
	centerG[1]=bl.position[1];

	auto tunel(-(generateDFCylinder(rTunel, makeAVec(l,0.,0.), center) &
	             generateDFPlane(makeAVec(0.,-1.,0.), centerG)));

	return normalize(tunel, dx);
}


int main()
{
	Param dx(0.5);
	Param dt(1.);
	Param nu(.01);

	
	AVec<int> size(makeAVec(150., 37., 53.2)*(1./dx.v()));
	asl::Block bl(size,dx.v(),makeAVec(-10.,-1.63, -31.6));
	
	Param nuNum(nu.v()*dt.v()/dx.v()/dx.v());
	
	std::cout<<"Flow: Data initialization...";


	auto object(asl::readSurf("train1.stl", bl));
	
	asl::Block block(object->getInternalBlock());

	auto tunelMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(tunelMap, generateTunel(block));

	auto forceField(asl::generateDataContainerACL_SP<FlT>(block, 3, 1u));
	asl::initData(forceField, makeAVec(0.,0.,0.));
	
	std::cout<<"Finished"<<endl;
	
	std::cout<<"Flow: Numerics initialization...";

	asl::SPLBGK lbgk(new asl::LBGK(block, 
				               acl::generateVEConstant(FlT(nuNum.v())),  
	        			       &asl::d3q15()));
	
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(.1,.0,.0));

	auto vfTunel(asl::generatePFConstant(makeAVec(0.1,0.,0.)));

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcV;

	bc.push_back(generateBCVelocity(lbgk, vfTunel, tunelMap));
	bcV.push_back(generateBCVelocityVel(lbgk, vfTunel, tunelMap));
//	bcV.push_back(generateBCNoSlipRho(lbgk, tunelMap));
	bc.push_back(generateBCNoSlip(lbgk,  object));
	bcV.push_back(generateBCNoSlipVel(lbgk, object));
//	bcV.push_back(generateBCNoSlipRho(lbgk, object));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1., makeAVec(0.1,0.,0.), {asl::X0, asl::XE}));

	initAll(bc);
	initAll(bcV);

	auto computeForce(generateComputeSurfaceForce(lbgk, forceField, object));
	computeForce->init();
	

	std::cout<<"Finished"<<endl;
	std::cout<<"Computing...";
	asl::Timer timer;

	asl::WriterVTKXML writer("flowAroundObjRes");
	writer.addScalars("map", *object);
	writer.addScalars("tunel", *tunelMap);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *lbgk->getVelocity());
	writer.addVector("force", *forceField);

	executeAll(bc);
	executeAll(bcV);
	computeForce->execute();

	writer.write();

	timer.start();
	for(unsigned int i(1); i < 2001; ++i)
	{
		lbgk->execute();
		executeAll(bc);
		if(!(i%200))
		{
			cout<<i<<endl;
			executeAll(bcV);
			computeForce->execute();
			writer.write();
		}
	}
	timer.stop();
	
	std::cout<<"Finished"<<endl;	

	cout << "time=" << timer.getTime() << "; clockTime="
		<< timer.getClockTime()	<< "; load=" 
		<< timer.getProcessorLoad() * 100 << "%" << endl;

	std::cout<<"Output...";
	std::cout<<"Finished"<<endl;	
	std::cout<<"Ok"<<endl;

	return 0;
}
