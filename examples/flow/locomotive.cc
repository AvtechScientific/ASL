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
	\example locomotive.cc
	Required input file: [locomotive.stl](http://asl.org.il/input_data/locomotive.stl)
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

	// Set length of the tunnel to the length (X size) of the block
	double l(bl.getBPosition()[0] - bl.position[0] + bl.dx);
	// Set radius of the tunnel to the ca. half of the block's height (Z size)
	double rTunnel((bl.getBPosition()[2] - bl.position[2]) / 2.1);

	// Center of the tunnel (described as cylinder cut by a plane)
	asl::AVec<> center(.5 * (bl.getBPosition() + bl.position));
	center[1] = bl.position[1] + .25 * rTunnel;

	// Center of the ground plane (that cuts the cylinder) 
	asl::AVec<> centerG(center); 
	centerG[1] = bl.position[1];

	/* DF = DistanceFunction (part of the geometrical module of ASL)
	 1. Genarate cylinder
	 2. Generate ground plane
	 3. Conjunction of the cylinder and the plane ('&' - operator)
	 4. Space inversion ('-' - operator) */
	auto tunnel(-(generateDFCylinder(rTunnel, makeAVec(l, 0., 0.), center) &
	             generateDFPlane(makeAVec(0., -1., 0.), centerG)));

	// Normalize DistanceFunction to the range [-1; 1]
	return normalize(tunnel, bl.dx);
}


int main(int argc, char* argv[])
{
	/* Convenience facility to manage simulation parameters (and also
	 hardware parameters defining platform and device for computations)
	 through command line and/or parameters file.
	 See `locomotive --help` for more information */
	asl::ApplicationParametersManager appParamsManager("locomotive",
	                                                   "1.0");

	/* Important: declare Parameters only after declaring
	 ApplicationParametersManager instance because each Parameter adds itself
	 to it automatically!
	 0.08 - default value; will be used if nothing else is provided during
	 runtime through command line or parameters file.
	 "dx" - option key; is used to specify this parameter through command line
	 and/or parameters file, like `locomotive --dx 0.05`
	 "space step" - option description; is used in the help output:
	 `locomotive -h` and as comment on parameters file generation:
	 `locomotive -g ./defaultParameters.ini`
	 "m" - parameter units; is used to complement the option description mentioned
	 above. Might be used for automatic unit conversion in future (to this end
	 it is recommended to use the notation of the Boost::Units library). */
	asl::Parameter<FlT> dx(0.08, "dx", "space step", "m");
	asl::Parameter<FlT> dt(1., "dt", "time step", "s");
	asl::Parameter<FlT> nu(.001, "nu", "kinematic viscosity", "m^2/s");
	asl::Parameter<unsigned int> iterations(10001, "iterations", "iterations number");
	asl::Parameter<string> input("input", "path to the geometry input file");

	/* Load previously declared Parameters from command line and/or
	parameters file. Use default values if neither is provided. */
	appParamsManager.load(argc, argv);

	/* Set the size of the block to 40x10x15 m (in accordance with the
	 locomotive size read later on from the input file) */ 
	AVec<int> size(makeAVec(40., 10., 15.) * (1. / dx.v()));

	/* Create block and shift it in accordance with the
	 position of the locomotive in the input file */
	asl::Block bl(size, dx.v(), makeAVec(-30., 8.58, 1.53));

	// Define dimensionless viscosity value
	FlT nuNum(nu.v() * dt.v() / dx.v() / dx.v());
	
	cout << "Data initialization... " << flush;

	// Read geometry of the locomotive from the file `locomotive.stl`
	auto locomotive(asl::readSurface(input.v(), bl));

	// Create block for further use
	asl::Block block(locomotive->getInternalBlock());

	// Generate memory data container for the tunnel
	auto tunnelMap(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	// Place generated geometry of the tunnel into the tunnel data container
	asl::initData(tunnelMap, generateTunnel(block));

	// Data container for air friction field
	auto forceField(asl::generateDataContainerACL_SP<FlT>(block, 3, 1u));
	// Initialization
	asl::initData(forceField, makeAVec(0., 0., 0.));
	
	cout << "Finished" << endl;

	cout << "Numerics initialization... " << flush;

	// NOTE: the problem is considered in the reference frame related to the locomotive

	// Generate numerical method for air flow - LBGK (lattice Bhatnagar–Gross–Krook)
	asl::SPLBGK lbgk(new asl::LBGKTurbulence(block, 
	                                         acl::generateVEConstant(FlT(nu.v())),  
	                                         &asl::d3q15()));
	lbgk->init();
	// Generate an instance for LBGK data initialization
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	// Initialize the LBGK internal data with the flow velocity of (0.1, 0, 0) in [lattice units]
	lbgkUtil->initF(acl::generateVEConstant(.1, .0, .0));

	auto vfTunnel(asl::generatePFConstant(makeAVec(0.1, 0., 0.)));

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcV;

	// Generate boundary conditions for the tunnel geometry. Constant velocity BC
	bc.push_back(generateBCVelocity(lbgk, vfTunnel, tunnelMap));
	// Generate boundary conditions for the tunnel geometry. Constant velocity BC
	// This BC is used for visualization.
	bcV.push_back(generateBCVelocityVel(lbgk, vfTunnel, tunnelMap));
	bcV.push_back(generateBCNoSlipRho(lbgk, tunnelMap));

	// Generate boundary conditions for the locomotive geometry. Non-slip BC
	bc.push_back(generateBCNoSlip(lbgk,  locomotive));
	bcV.push_back(generateBCNoSlipVel(lbgk, locomotive));
	bcV.push_back(generateBCNoSlipRho(lbgk, locomotive));

	// Generate constant presure BC for in and out planes of the tunnel
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1.,
	                                                makeAVec(0.1, 0., 0.),
	                                                {asl::X0, asl::XE}));

	// Initialization and building of all BC
	initAll(bc);
	initAll(bcV);

	// Generate a numerical method for computation of the air force field that acts on the locomotive
	auto computeForce(generateComputeSurfaceForce(lbgk, forceField, locomotive));
	computeForce->init();
	
	cout << "Finished" << endl;
	cout << "Computing..." << endl;
	asl::Timer timer;

	// Initialization of the output system
	// Write the output to the directory containing the input parameters file (default "./")
	asl::WriterVTKXML writer(appParamsManager.getDir() + "locomotive");
	writer.addScalars("map", *locomotive);
	writer.addScalars("tunnel", *tunnelMap);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *lbgk->getVelocity());
	writer.addVector("force", *forceField);

	// Execute all BC
	executeAll(bc);
	executeAll(bcV);
	computeForce->execute();

	// First data output
	writer.write();

	timer.start();
	// Iteration loop 
	for (unsigned int i(1); i < iterations.v(); ++i)
	{
		// One iteration (timestep) of bulk numerical procedure
		lbgk->execute();
		// Execution of the BC procedures
		executeAll(bc);
		// Output and analysis scope
		if (!(i%1000))
		{
			cout << i << endl;
			// Execution of the visualization BC procedures
			executeAll(bcV);
			// Computation of the force field
			computeForce->execute();
			// Data writing
			writer.write();
		}
	}
	timer.stop();
	
	cout << "Finished" << endl;	

	cout << "Computation statistic:" << endl;
	cout << "time = " << timer.getTime() << "; clockTime = "
		 << timer.getClockTime() << "; load = "
		 << timer.getProcessorLoad() * 100 << "%" << endl;

	return 0;
}