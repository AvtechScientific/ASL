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
	\example acousticWaves.cc
 */

#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslFDElasticity.h>
#include <num/aslFDElasticityBC.h>
#include <num/aslBasicBC.h>
#include <utilities/aslTimer.h>
#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>


typedef float FlT;
//typedef asl::UValue<FlT> Param;

class Parameters
{
  private:
		void init();

  public:
		asl::ApplicationParametersManager appParamsManager;

		asl::Block::DV size;

		asl::Parameter<double> dx;
		asl::Parameter<double> bulkModulus;
		asl::Parameter<double> shearModulus;
		asl::Parameter<double> rho;
		asl::Parameter<double> tubeL;
		asl::Parameter<double> tubeDEx;
		asl::Parameter<double> tubeDIn;
		asl::Parameter<double> hole1Pos;
		asl::Parameter<double> hole2Pos;
		asl::Parameter<double> hole1D;
		asl::Parameter<double> hole2D;

		asl::Parameter<double> tSimulation;
		asl::Parameter<double> tOutput;

		asl::UValue<double> dt;
		asl::UValue<double> bulkMNum;
		asl::UValue<double> shearMNum;
		
		void load(int argc, char * argv[]);
		Parameters();
		void updateNumValues();
};


Parameters::Parameters():
	appParamsManager("acousticWaves", "0.1"),
	size(3),
	dx(1e-3,"dx", "dx"),
	bulkModulus(160e9,"bulk_modulus", "bulk modulus"),
	shearModulus(79e9,"shear_modulus", "shear modulus"),
	rho(7800,"rho", "density"),
	tubeL(.2,"tube_length", "pipe length" "m"),
	tubeDEx(0.021, "tube_diameter_external", "external pipe diameter" "m"),
//	tubeDIn(0.0157,"tube_diameter_internal", "internal pipe diameter" "m"),
	tubeDIn(0.0107,"tube_diameter_internal", "internal pipe diameter" "m"),
	hole1Pos(0.1,"hole_1_position", "position of first hole" "m"),
	hole2Pos(0.15,"hole_2_position", "position of second hole" "m"),
	hole1D(15e-3,"hole_1_diameter", "diameter of first hole" "m"),
	hole2D(15e-3,"hole_2_diameter", "diameter of second hole" "m"),
	tSimulation(8e-5, "simulation_time", "simulation time"),
	tOutput(1e-6, "output_interval", "output interval")
{
}


void Parameters::load(int argc, char * argv[])
{
	appParamsManager.load(argc, argv);

	init();
}


void Parameters::updateNumValues()
{
	double vs(sqrt((bulkModulus.v()+2.*shearModulus.v())/rho.v()));
	dt=dx.v()/vs*.1;
	cout << vs << "; " << dx.v() << "; " << dt.v() << endl;
	bulkMNum = bulkModulus.v()/rho.v()/dx.v()/dx.v();
	shearMNum = shearModulus.v()/rho.v()/dx.v()/dx.v();
	size = asl::makeAVec(tubeL.v() / dx.v() + 1, 
	                     tubeDEx.v() / dx.v() + 1, 
	                     tubeDEx.v() / dx.v() + 1);
}


void Parameters::init()
{
//	if (tubeD.v() < pumpD.v())
//		asl::errorMessage("Tube's diameter is smaller than pump's diameter");
	updateNumValues();
}

asl::SPDistanceFunction generatePipe(asl::Block & block, Parameters &params)
{
	asl::SPDistanceFunction pipeGeometry;
	asl::AVec<double> orientation(asl::makeAVec(1., 0., 0.));
	asl::AVec<double> lVec(asl::makeAVec(params.tubeL.v()+2.*params.dx.v(), 0., 0.));
	asl::AVec<double> h1Orientation(asl::makeAVec(0., 1., 0.));
	asl::AVec<double> h2Orientation(asl::makeAVec(0., 0., 1.));
	asl::AVec<double> center(asl::AVec<double>(params.size)*.5*params.dx.v());
	double wallMid((params.tubeDEx.v()+params.tubeDIn.v())*.25);
	double wallTh((params.tubeDEx.v())*.5);
	asl::AVec<double> h1Center(center - (center*orientation)*orientation + 
	                           params.hole1Pos.v()*orientation + 
	                           h1Orientation*wallMid);
	asl::AVec<double> h2Center(center - (center*orientation)*orientation + 
	                           params.hole2Pos.v()*orientation + 
	                           h2Orientation*wallMid);

	pipeGeometry = asl::generateDFCylinder(params.tubeDEx.v() / 2., lVec, center) &
					(-asl::generateDFCylinderInf(params.tubeDIn.v() / 2., orientation, center));
	pipeGeometry = pipeGeometry &  
		           (-asl::generateDFCylinder(params.hole1D.v() / 2., h1Orientation * wallTh, h1Center));
	pipeGeometry = pipeGeometry &  
		           (-asl::generateDFCylinder(params.hole2D.v() / 2., h2Orientation * wallTh, h2Center));	
	return asl::normalize(-pipeGeometry, params.dx.v());
}

asl::AVec<float> getAmplitude(double it)
{
	float a(it<200. ? 1.-cos(it*6.28/200.) : 0);
	return asl::makeAVec(a,0.f,0.f);
}


int main(int argc, char* argv[])
{
	Parameters params;
	params.load(argc, argv);
		
	std::cout << "Data initialization... " << flush;

	asl::Block block(params.size, params.dx.v());
	auto displacement(asl::generateDataContainerACL_SP<FlT>(block, 3, 1u));
	asl::initData(displacement, asl::makeAVec(0.,0.,0.));

	auto mapMem(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(mapMem, generatePipe(block, params));

	
	asl::WriterVTKXML writer(params.appParamsManager.getDir() + "acousticWaves");
	writer.addScalars("map", *mapMem);
	writer.addVector("displacement", *displacement);
	writer.write();

	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... " << flush;

	auto elasticity(generateFDElasticityRelax(displacement,
	                                          params.bulkMNum.v(),
	                                          params.shearMNum.v(), 
	                                          params.dt.v(),
	                                          &asl::d3q19()));
/*	auto elasticity(generateFDElasticity(displacement,
	                                     params.bulkMNum.v(),
	                                     params.shearMNum.v(), 
	                                     params.dt.v(),
	                                     &asl::d3q19()));*/
	elasticity->setDumpingFactor(acl::generateVEConstant(.9999));
	elasticity->init();


	std::vector<asl::SPNumMethod> bc;

	bc.push_back(generateBCZeroStress(elasticity, mapMem));
//    bc.push_back(generateBCConstantGradient(displacement,asl::makeAVec(0.,0.,0.),mapMem,&asl::d3q19()));
//    bc.push_back(generateBCConstantGradient(elasticity->getPressureData(),0.,mapMem,&asl::d3q19()));
//	bc.push_back(generateBCZeroStressP(elasticity, mapMem));
	asl::UValue<asl::AVec<float>> pres(asl::makeAVec(0.f,0.f,0.f));
	bc.push_back(asl::generateBCConstantValue(displacement, pres, {asl::X0}));

	initAll(bc);

	std::cout << "Finished" << endl;
	std::cout << "Computing..." << endl;
	asl::Timer timer;

	executeAll(bc);

	timer.start();
	double tOutPrev(0);
	cout << params.dt.v() << endl;
	for (double t(0); t < params.tSimulation.v(); t+=params.dt.v())
	{
		elasticity->execute();
		pres=getAmplitude(t/params.dt.v());
		executeAll(bc);
		if(t - params.tOutput.v()>=tOutPrev)
		{
			timer.stop();
			tOutPrev=t;
			cout << t << "/" << params.tSimulation.v() << "; expected left time: " <<  
				  timer.getLeftTime(t/params.tSimulation.v())  << endl;
			writer.write();
			timer.resume();

		}
	}
	timer.stop();
	
	std::cout << "Finished" << endl;	

	cout << "time=" << timer.getTime() << "; clockTime="
		 <<  timer.getClockTime()	 <<  "; load=" 
		 <<  timer.getProcessorLoad() * 100 << "%" << endl;

	std::cout << "Output...";
	std::cout << "Finished" << endl;	
	std::cout << "Ok" << endl;

	return 0;
}
