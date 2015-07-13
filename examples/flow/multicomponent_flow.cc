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
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLBGK.h>
#include <num/aslLBGKBC.h>
#include <utilities/aslTimer.h>
#include <num/aslFDAdvectionDiffusion.h>
#include <num/aslBasicBC.h>

// typedef to switch to double precision
//typedef double FlT;
typedef float FlT;

using asl::AVec;
using asl::makeAVec;

class Parameters
{
  private:
		void init();

  public:
		asl::ApplicationParametersManager appParamsManager;

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

		asl::Parameter<double> component2InVel;
		asl::Parameter<double> component1InVel;
		asl::Parameter<double> component3InVel;
		
		
		void load(int argc, char * argv[]);
		Parameters();
		void updateNumValues();
};


Parameters::Parameters():
	appParamsManager("multicomponent_flow", "0.1"),
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
	component1InVel(0.16, "component1_in_velocity", "flow velocity in the component1 input"),
	component2InVel(0.08, "component2_in_velocity", "flow velocity in the component2 input"),
	component3InVel(0.1, "component3_in_velocity", "flow velocity in the component3 input")
{
}


void Parameters::load(int argc, char * argv[])
{
	appParamsManager.load(argc, argv);
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

// Generate geometry of the mixer (cross-coupled pipes)
asl::SPDistanceFunction generateMixer(asl::Block & block, Parameters &params)
{
	asl::SPDistanceFunction mixerGeometry;
	asl::AVec<double> orientation(asl::makeAVec(0., 0., 1.));
	asl::AVec<double> center(asl::AVec<double>(params.size) * .5 * params.dx.v());

	mixerGeometry = generateDFCylinderInf(params.tubeD.v() / 2., orientation,
	                                      center);

	orientation[1] = 1.0;
	orientation[2] = 0.0;
	center[2] = params.pumpD.v() * 1.5;
	mixerGeometry = mixerGeometry | generateDFCylinderInf(params.pumpD.v() / 2.,
	                                                      orientation, center);

	return asl::normalize(-(mixerGeometry) | asl::generateDFInBlock(block, 0),
	                      params.dx.v());
}

int main(int argc, char *argv[])
{
	Parameters params;
	params.load(argc, argv);
	
	cout << "Data initialization..." << endl;

	asl::Block block(params.size, params.dx.v());

	auto mcfMapMem(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(mcfMapMem, generateMixer(block, params));

	auto component1Frac(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(component1Frac, 0);
	auto component3Frac(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(component3Frac, 0);
	
	
	cout << "Finished" << endl;

	cout << "Numerics initialization..." << endl;

	auto templ(&asl::d3q15());	
	
	asl::SPLBGK lbgk(new asl::LBGK(block,
	                               acl::generateVEConstant(FlT(params.nuNum.v())),
	                               templ));
	
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(.0, .0, .0));

	auto flowVel(lbgk->getVelocity());
	auto nmcomponent1(asl::generateFDAdvectionDiffusion(component1Frac, 0.01,
	                                                    flowVel, templ, true));
	nmcomponent1->init();
	auto nmcomponent3(asl::generateFDAdvectionDiffusion(component3Frac, 0.01,
	                                                    flowVel, templ));
	nmcomponent3->init();

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcV;
	std::vector<asl::SPNumMethod> bcDif;
	
	bc.push_back(generateBCNoSlip(lbgk, mcfMapMem));
	bc.push_back(generateBCConstantPressure(lbgk, 1., {asl::ZE}));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1., 
	                                                makeAVec(0., 0., params.component2InVel.v()),
	                                                {asl::Z0}));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1., 
	                                                makeAVec(0., -params.component1InVel.v(), 0.),
	                                                {asl::YE}));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1., 
	                                                makeAVec(0., params.component3InVel.v(), 0.),
	                                                {asl::Y0}));
	
	bcDif.push_back(generateBCNoSlipVel(lbgk, mcfMapMem));
	bc.push_back(generateBCConstantGradient(component1Frac, 0., mcfMapMem, templ));
	bc.push_back(generateBCConstantGradient(component3Frac, 0., mcfMapMem, templ));
	bc.push_back(generateBCConstantValue(component1Frac, 1., {asl::YE}));
	bc.push_back(generateBCConstantValue(component3Frac, 0., {asl::YE, asl::Z0, asl::ZE}));
	bc.push_back(generateBCConstantValue(component1Frac, 0., {asl::Y0, asl::Z0, asl::ZE}));
	bc.push_back(generateBCConstantValue(component3Frac, 1., {asl::Y0}));
//	bc.push_back(generateBCConstantGradient(component1Frac, 0.,templ, {asl::ZE}));
//	bc.push_back(generateBCConstantGradient(component3Frac, 0.,templ, {asl::ZE}));
	
	initAll(bc);
	initAll(bcDif);
	initAll(bcV);

	cout << "Finished" << endl;
	cout << "Computing..." << endl;
	asl::Timer timer;

	asl::WriterVTKXML writer("multicomponent_flow");
	writer.addScalars("map", *mcfMapMem);
	writer.addScalars("component1", *component1Frac);
	writer.addScalars("component3", *component3Frac);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *flowVel);

	executeAll(bc);
	executeAll(bcDif);
	executeAll(bcV);

	writer.write();

	timer.start();
	for (unsigned int i(1); i < 10001; ++i)
	{
		lbgk->execute();
		executeAll(bcDif);
		nmcomponent1->execute();
		nmcomponent3->execute();
		executeAll(bc);
		
		if (!(i%100))
		{
			timer.stop();
			cout << i << "/10000; time left (expected): " <<  timer.getLeftTime(double(i)/10000.)  << endl;
			executeAll(bcV);
			writer.write();
			timer.resume();
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
