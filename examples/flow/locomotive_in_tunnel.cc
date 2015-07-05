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
	\example locomotive_in_tunnel.cc
 */

#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <math/aslPositionFunction.h>
#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLBGK.h>
#include <num/aslLBGKBC.h>
#include <utilities/aslTimer.h>
#include <readers/aslVTKFormatReaders.h>


// typedef to switch to double precision
//typedef double FlT;
typedef float FlT;

using asl::AVec;
using asl::makeAVec;

// Generate geometry of the tunnel
asl::SPDistanceFunction generateTunnel(asl::Block & bl)
{

	double l(bl.getBPosition()[0] - bl.position[0] + bl.dx);
	double rTunnel((bl.getBPosition()[2] - bl.position[2]) / 2.1);

	double dx(bl.dx);

	asl::AVec<int> size(bl.getSize());

	asl::AVec<> center(.5 * (bl.getBPosition() + bl.position)); 
	center[1] = bl.position[1] + .25 * rTunnel;
	asl::AVec<> centerG(center); 
	centerG[1] = bl.position[1];

	auto tunnel(-(generateDFCylinder(rTunnel, makeAVec(l, 0., 0.), center) &
	             generateDFPlane(makeAVec(0., -1., 0.), centerG)));

	return normalize(tunnel, dx);
}


int main(int argc, char* argv[])
{
	/* Convenience facility to manage simulation parameters (and also
	hardware parameters - platform/device to run the application on)
	through command line and/or parameters file.
	See `locomotive_in_tunnel -h` for more information */
	asl::ApplicationParametersManager appParamsManager("locomotive_in_tunnel",
	                                                   "1.0");

	/* Important: declare Parameters only after declaring
	ApplicationParametersManager instance because each Parameter adds itself
	to it automatically */
	asl::Parameter<FlT> dx(0.08, "dx", "space step", "m");
	asl::Parameter<FlT> dt(1., "dt", "time step", "s");
	asl::Parameter<FlT> nu(.001, "nu", "viscosity", "Pa*s");

	/* Load previously declared Parameters from command line and/or
	parameters file. Use default values if neither is provided. */
	appParamsManager.load(argc, argv);

	AVec<int> size(makeAVec(40., 10., 15.) * (1. / dx.v()));
	asl::Block bl(size, dx.v(), makeAVec(-30., 8.58, 1.53));
	
	asl::UValue<FlT> nuNum(nu.v() * dt.v() / dx.v() / dx.v());
	
	std::cout << "Data initialization... ";

	auto locomotive(asl::readSurface("locomotive.stl", bl));
	
	asl::Block block(locomotive->getInternalBlock());

	auto tunnelMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(tunnelMap, generateTunnel(block));

	auto forceField(asl::generateDataContainerACL_SP<FlT>(block, 3, 1u));
	asl::initData(forceField, makeAVec(0., 0., 0.));
	
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... ";

	asl::SPLBGK lbgk(new asl::LBGKTurbulence(block, 
	                                         acl::generateVEConstant(FlT(nu.v())),  
	                                         &asl::d3q15()));
	
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(.1, .0, .0));

	auto vfTunnel(asl::generatePFConstant(makeAVec(0.1, 0., 0.)));

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcV;

	bc.push_back(generateBCVelocity(lbgk, vfTunnel, tunnelMap));
	bcV.push_back(generateBCVelocityVel(lbgk, vfTunnel, tunnelMap));
//	bcV.push_back(generateBCNoSlipRho(lbgk, tunnelMap));
	bc.push_back(generateBCNoSlip(lbgk,  locomotive));
	bcV.push_back(generateBCNoSlipVel(lbgk, locomotive));
//	bcV.push_back(generateBCNoSlipRho(lbgk, locomotive));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1.,
	                                                makeAVec(0.1, 0., 0.),
	                                                {asl::X0, asl::XE}));

	initAll(bc);
	initAll(bcV);

	auto computeForce(generateComputeSurfaceForce(lbgk, forceField, locomotive));
	computeForce->init();
	

	std::cout << "Finished" << endl;
	std::cout << "Computing... ";
	asl::Timer timer;

	asl::WriterVTKXML writer(appParamsManager.getDir() + "locomotive_in_tunnel");
	writer.addScalars("map", *locomotive);
	writer.addScalars("tunnel", *tunnelMap);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *lbgk->getVelocity());
	writer.addVector("force", *forceField);

	executeAll(bc);
	executeAll(bcV);
	computeForce->execute();

	writer.write();

	timer.start();
	for (unsigned int i(1); i < 20001; ++i)
	{
		lbgk->execute();
		executeAll(bc);
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
		 << timer.getClockTime() << "; load=" 
		 << timer.getProcessorLoad() * 100 << "%" << endl;

	std::cout << "Output... ";
	std::cout << "Finished" << endl;	
	std::cout << "Ok" << endl;

	return 0;
}