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
	\example multicomponent_flow.cc 	 
 */

#include <utilities/aslParametersManager.h>
#include <math/aslVectors.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <data/aslDataWithGhostNodes.h>
#include <aslGenerators.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLBGK.h>
#include <num/aslLBGKBC.h>
#include "utilities/aslTimer.h"
#include "acl/aclUtilities.h"
#include <num/aslFDAdvectionDiffusion.h>
#include <num/aslBasicBC.h>


typedef float FlT;
//typedef double FlT;
typedef asl::UValue<double> Param;

using asl::AVec;
using asl::makeAVec;

class Parameters
{
  private:
		void init();

  public:
		asl::ParametersManager parametersManager;
		string folder;

		asl::Block::DV size;

		asl::Parameter<double> dx;
		asl::Parameter<double> dt;

		asl::Parameter<double> tSimulation;
		asl::Parameter<double> tOutput;

		asl::Parameter<double> nu;
		asl::UValue<double> nuNum;
		
		asl::Parameter<double> tubeL;
		asl::Parameter<double> tubeD;
		asl::Parameter<double> pumpL;
		asl::Parameter<double> pumpD;

		asl::Parameter<double> oilInVel;
		asl::Parameter<double> waterInVel;
		asl::Parameter<double> gasInVel;
		
		
		void load(int argc, char * argv[],
		          string programName,
		          string programVersion);
		Parameters();
		void updateNumValues();
};


Parameters::Parameters():
	size(3),
	dx(0.0005, "dx", "space step"),
	dt(1., "dt", "time step"),
	tSimulation(2e-3, "simulation_time", "simulation time"),
	tOutput(1e-4, "output_interval", "output interval"),
	nu(4e-8/1.6, "nu", "viscosity"),
	tubeL(0.25, "tubeL", "tube's length"),
	tubeD(0.05, "tubeD", "tube's diameter"),
	pumpL(0.025, "pumpL", "pump's length"),
	pumpD(0.03, "pumpD", "pump's diameter"),
	oilInVel(0.08, "oil_in_velocity", "flow velocity in the oil input"),
	waterInVel(0.16, "water_in_velocity", "flow velocity in the water input"),
	gasInVel(0.1, "gas_in_velocity", "flow velocity in the gas input")
{
}


void Parameters::load(int argc, char * argv[],
                      string programName,
                      string programVersion)
{
	parametersManager.load(argc, argv, programName, programVersion);
	folder = parametersManager.getFolderWithSlash();

	init();
}


void Parameters::updateNumValues()
{
	nuNum = nu.v() * dt.v() / dx.v() / dx.v();
	size[0] = tubeD.v() / dx.v() + 1;
	size[1] = (tubeD.v() + 2 * pumpL.v()) / dx.v() + 1;
	size[2] = tubeL.v() / dx.v() + 1;
}


void Parameters::init()
{
	if (tubeD.v() < pumpD.v())
		asl::errorMessage("Tube's diameter is smaller than pump's diameter");

	updateNumValues();
}

// generate geometry of the mixer
asl::SPDistanceFunction generateMixer(asl::Block & block, Parameters &params)
{
	asl::SPDistanceFunction mixerGeometry;
	asl::AVec<double> orientation(asl::makeAVec(0., 0., 1.));
	asl::AVec<double> center(asl::AVec<double>(params.size)*.5*params.dx.v());

	mixerGeometry = generateDFCylinderInf(params.tubeD.v() / 2., orientation, center);

	orientation[1] = 1.0;
	orientation[2] = 0.0;
	center[2]=params.pumpD.v() * 1.5;
	mixerGeometry = mixerGeometry | generateDFCylinderInf(params.pumpD.v() / 2., orientation, center);

	return asl::normalize(-(mixerGeometry) | asl::generateDFInBlock(block, 0), params.dx.v());
}

int main(int argc, char *argv[])
{
	Parameters params;
	params.load(argc, argv, "multicomponent_flow", "0.1");
	
	std::cout<<"Flow: Data initialization...";

	asl::Block block(params.size, params.dx.v());

	auto mcfMapMem(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(mcfMapMem, generateMixer(block, params));

	auto waterFrac(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(waterFrac, 0);
	auto gasFrac(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(gasFrac, 0);
	
	
	std::cout<<"Finished"<<endl;
	
	std::cout<<"Flow: Numerics initialization...";

	auto templ(&asl::d3q15());	
	
	asl::SPLBGK lbgk(new asl::LBGK(block,
	                               acl::generateVEConstant(FlT(params.nuNum.v())),
	                               templ));
	
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(.0, .0, .0));

	auto flowVel(lbgk->getVelocity());
	auto nmWater(asl::generateFDAdvectionDiffusion(waterFrac, 0.01, flowVel, templ, true));
	nmWater->init();
	auto nmGas(asl::generateFDAdvectionDiffusion(gasFrac, 0.01, flowVel, templ));
	nmGas->init();

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcV;
	std::vector<asl::SPNumMethod> bcDif;
	
	bc.push_back(generateBCNoSlip(lbgk, mcfMapMem));
	bc.push_back(generateBCConstantPressure(lbgk,1.,{asl::ZE}));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1., 
	                                                makeAVec(0.,0.,params.oilInVel.v()),
	                                                {asl::Z0}));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1., 
	                                                makeAVec(0.,-params.waterInVel.v(),0.),
	                                                {asl::YE}));
	bc.push_back(generateBCConstantPressureVelocity(lbgk,1., 
	                                                makeAVec(0.,params.gasInVel.v(),0.),
	                                                {asl::Y0}));
	
	bcDif.push_back(generateBCNoSlipVel(lbgk, mcfMapMem));
	bc.push_back(generateBCConstantGradient(waterFrac, 0., mcfMapMem, templ));
	bc.push_back(generateBCConstantGradient(gasFrac, 0., mcfMapMem, templ));
	bc.push_back(generateBCConstantValue(waterFrac, 1., {asl::YE}));
	bc.push_back(generateBCConstantValue(gasFrac, 0., {asl::YE, asl::Z0, asl::ZE}));
	bc.push_back(generateBCConstantValue(waterFrac, 0., {asl::Y0, asl::Z0, asl::ZE}));
	bc.push_back(generateBCConstantValue(gasFrac, 1., {asl::Y0}));
//	bc.push_back(generateBCConstantGradient(waterFrac, 0.,templ, {asl::ZE}));
//	bc.push_back(generateBCConstantGradient(gasFrac, 0.,templ, {asl::ZE}));
	
	initAll(bc);
	initAll(bcDif);
	initAll(bcV);

	std::cout<<"Finished"<<endl;
	std::cout<<"Computing..."<<endl;
	asl::Timer timer;

	asl::WriterVTKXML writer("multicomponent_flow");
	writer.addScalars("map", *mcfMapMem);
	writer.addScalars("water", *waterFrac);
	writer.addScalars("gas", *gasFrac);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *flowVel);

	executeAll(bc);
	executeAll(bcDif);
	executeAll(bcV);

	writer.write();

	timer.start();
	for(unsigned int i(1); i < 10001; ++i)
	{
		lbgk->execute();
		executeAll(bcDif);
		nmWater->execute();
		nmGas->execute();
		executeAll(bc);
		
		if(!(i%100))
		{
			timer.stop();
			cout<<i<<"/10000; expected left time: "<< timer.getLeftTime(double(i)/10000.) <<endl;
			executeAll(bcV);
			writer.write();
			timer.resume();
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
