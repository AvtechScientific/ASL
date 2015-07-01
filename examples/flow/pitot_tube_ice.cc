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
	\example pitot_tube_ice.cc
    Icing process in the Pitot tube
 */

#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <aslGenerators.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLBGK.h>
#include <num/aslLBGKBC.h>
#include <utilities/aslTimer.h>
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
		asl::ApplicationParametersManager appParamsManager;

		asl::Block::DV size;

		asl::Parameter<double> dx;
		asl::Parameter<double> dt;

		asl::Parameter<double> tSimulation;
		asl::Parameter<double> tOutput;

		asl::Parameter<double> nu;
		asl::UValue<double> nuNum;
		
		asl::Parameter<double> rIn;
		asl::Parameter<double> rEx;
		asl::Parameter<double> lCyl;
		asl::Parameter<double> lCone;

		asl::Parameter<double> temperature;
		asl::Parameter<double> humidity;
		asl::Parameter<double> flowVel;
		
		void load(int argc, char * argv[]);
		Parameters();
		void updateNumValues();
};


Parameters::Parameters():
	appParamsManager("pitot_tube_ice", "0.1"),
	size(3),
	dx(0.000125, "dx", "space step"),
	dt(1., "dt", "time step"),
	tSimulation(2e-3, "simulation_time", "simulation time"),
	tOutput(1e-4, "output_interval", "output interval"),
	nu(6.25e-10/4., "nu", "viscosity"),
	rIn(0.0015, "r_in", "Internal radius, m"),
	rEx(0.005, "r_ex", "External radius, m"),
	lCyl(0.002, "l_cyl", "Length of cylindric part, m"),
	lCone(0.02, "l_cone", "Length of conic part, m"),
	temperature(253, "temperature", "temperature, K"),
	humidity(.5, "humidity", "relative humidity, K"),
	flowVel(0.08, "flow_vel", "flow velocity")
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
	size[0] = 1.0*(lCyl.v() + lCone.v()) / dx.v() + 1;
	size[1] = rEx.v() * 2.5 / dx.v() + 1;
	size[2] = rEx.v() * 2.5 / dx.v() + 1;
}


void Parameters::init()
{
	if (rEx.v() < rIn.v())
		asl::errorMessage("External tube's diameter is smaller than internal one");

	updateNumValues();
}

// generate geometry
asl::SPDistanceFunction generateGeometry(asl::Block & block, Parameters &params)
{
	asl::SPDistanceFunction tubeGeometry;
	asl::AVec<double> orientation(asl::makeAVec(1., 0., 0.));
	asl::AVec<double> center(asl::AVec<double>(params.size)*.5*params.dx.v());
	auto centerCyl(center);
	centerCyl[0] = params.lCyl.v()*.45;
	auto centerHole(centerCyl+(params.lCone.v()*.6)*orientation);
	auto lHole(params.lCyl.v()+params.lCone.v());
		
	auto apexCone(centerCyl+orientation*(params.lCyl.v()*.49+params.lCone.v()));

	tubeGeometry = ((generateDFCylinder(params.rEx.v(), orientation*params.lCyl.v(), centerCyl) |
	                generateDFCone(params.rEx.v()*.98, -orientation*params.lCone.v(), apexCone)) &
		           -generateDFCylinder(params.rIn.v(), orientation*lHole, centerHole)) &
					generateDFPlane(orientation, apexCone-orientation*params.lCone.v()*.5);

	return asl::normalize(tubeGeometry, params.dx.v());
}

int main(int argc, char *argv[])
{
	Parameters params;
	params.load(argc, argv);
	
	std::cout << "Data initialization...";

	asl::Block block(params.size, params.dx.v());

	auto mcfMapMem(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(mcfMapMem, generateGeometry(block, params));

//	auto waterFrac(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
//	asl::initData(waterFrac, 0);
	
	
	std::cout << "Finished" << endl;
	
	std::cout << "Flow: Numerics initialization...";

	auto templ(&asl::d3q15());	
	
	asl::SPLBGK lbgk(new asl::LBGK(block,
	                               acl::generateVEConstant(FlT(params.nuNum.v())),
	                               templ));
	
	lbgk->init();
	asl::SPLBGKUtilities lbgkUtil(new asl::LBGKUtilities(lbgk));
	lbgkUtil->initF(acl::generateVEConstant(-0.9*params.flowVel.v(), params.flowVel.v()*.4, .0));

	auto flowVel(lbgk->getVelocity());
//	auto nmWater(asl::generateFDAdvectionDiffusion(waterFrac, 0.01, flowVel, templ, false));
//	nmWater->init();

	std::vector<asl::SPNumMethod> bc;
	std::vector<asl::SPNumMethod> bcV;
	std::vector<asl::SPNumMethod> bcDif;
	
	bc.push_back(generateBCNoSlip(lbgk, mcfMapMem));
	bc.push_back(generateBCConstantPressure(lbgk,1.,{asl::ZE}));
	bc.push_back(generateBCConstantPressureVelocity(lbgk, 1., 
	                                                makeAVec(-params.flowVel.v()*.9,params.flowVel.v()*.3,0.),
	                                                {asl::X0,asl::XE,asl::Y0,asl::YE,asl::Z0,asl::ZE}));
	
	bcDif.push_back(generateBCNoSlipVel(lbgk, mcfMapMem));
	bcV.push_back(generateBCNoSlipRho(lbgk, mcfMapMem));
//	bc.push_back(generateBCConstantGradient(waterFrac, 0., mcfMapMem, templ));
//	bc.push_back(generateBCConstantValue(waterFrac, 1., {asl::YE}));
//	bc.push_back(generateBCConstantValue(waterFrac, 0., {asl::Y0, asl::Z0, asl::ZE}));
	
	initAll(bc);
	initAll(bcDif);
	initAll(bcV);

	std::cout << "Finished" << endl;
	std::cout << "Computing..." << endl;
	asl::Timer timer;

	asl::WriterVTKXML writer(params.appParamsManager.getDir() + "pitot_tube");
	writer.addScalars("map", *mcfMapMem);
//	writer.addScalars("water", *waterFrac);
	writer.addScalars("rho", *lbgk->getRho());
	writer.addVector("v", *flowVel);

	executeAll(bc);
	executeAll(bcDif);
	executeAll(bcV);

	writer.write();

	timer.start();
	for(unsigned int i(1); i < 8001; ++i)
	{
		lbgk->execute();
		executeAll(bcDif);
//		nmWater->execute();
		executeAll(bc);
		
		if(!(i%800))
		{
			timer.stop();
			cout << i << "/8000; expected left time: " <<  timer.getLeftTime(double(i)/8000.)  << endl;
			executeAll(bcV);
			writer.write();
			timer.resume();
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
