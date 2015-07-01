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
	\example cubePoroelasticGravity.cc
 */

#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslFDPoroElasticity.h>
#include <num/aslFDElasticityBC.h>
#include <num/aslFDPoroElasticityBC.h>
#include <utilities/aslTimer.h>
#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <aslGeomInc.h>
#include <acl/aclUtilities.h>


typedef float FlT;
//typedef double FlT;
typedef asl::UValue<FlT> Param;


int main(int argc, char* argv[])
{
	asl::ApplicationParametersManager appParamsManager("cubePoroelasticGravity",
	                                                   "1.0");
	asl::Parameter<asl::AVec<int> > size("size", "size 3D");
	asl::Parameter<cl_float> dx("dx", "dx");
	asl::Parameter<cl_float> dt("dt", "dt");
	asl::Parameter<cl_float> bulkModulus("bulk_modulus", "bulk modulus");
	asl::Parameter<cl_float> shearModulus("shear_modulus", "shear modulus");
	asl::Parameter<cl_float> hydraulicConductivity("hydraulic_conductivity", "hydraulic conductivity");
	asl::Parameter<cl_float> rho("rho", "density");
	asl::Parameter<asl::AVec<FlT> > g("g", "gravity vector");

	asl::Parameter<unsigned int> tsim("num_iterations", "number of iterations");
	asl::Parameter<unsigned int> tout("num_it_out", "number of iterations between outputs");
	
	appParamsManager.load(argc, argv);

	Param bulkModulusNum(bulkModulus.v()/rho.v()/dx.v()/dx.v()*dt.v()*dt.v());
	Param shearModulusNum(shearModulus.v()/rho.v()/dx.v()/dx.v()*dt.v()*dt.v());

	asl::AVec<FlT> gNum(g.v()*dt.v()*dt.v());
		
	std::cout << "Data initialization... " << flush;

	asl::Block block(size.v(), dx.v());
	auto displacement(asl::generateDataContainerACL_SP<FlT>(block, 3, 1u));
	asl::initData(displacement, asl::makeAVec(0.,0.,0.));
	auto pressureL(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	asl::initData(pressureL, 0.);

	auto mapDF(asl::normalize(asl::generateDFInBlock(block, 0), dx.v()));
	auto map(asl::generateDataContainer_SP(block, mapDF, 1u, acl::typeToTypeID<FlT>()));
	auto mapX(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	initData(mapX, mapDF);

	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... " << flush;

	auto elasticity(generateFDPoroElasticity(displacement,
	                                         pressureL,
	                                         bulkModulusNum.v(),
	                                         shearModulusNum.v(),
	                                         hydraulicConductivity.v(),
	                                         &asl::d3q15()));
	elasticity->setForce(acl::generateVEConstant(gNum));
	elasticity->init();

	std::vector<asl::SPNumMethod> bc;
	asl::addBCZeroStress(bc, elasticity, mapX);
	asl::addBCRigidWall(bc, elasticity, {asl::X0});
	asl::initAll(bc);

	asl::WriterVTKXML writer(appParamsManager.getDir() + "cubePoroelasticGravity");
	writer.addScalars("map", *mapX);
	writer.addVector("displacement", *displacement);
	writer.addScalars("pressure", *elasticity->getPressureData());
	writer.addScalars("pressureL", *elasticity->getLiquidPressureData());
	
	std::cout << "Finished" << endl;
	std::cout << "Computing..." << endl;
	asl::Timer timer, timerBulk, timerBC;

	executeAll(bc);
	writer.write();

	timer.start();
	timerBulk.reset();
	timerBC.reset();
	for (unsigned int i(0); i < tsim.v(); ++i)
	{
		timerBulk.resume();
		elasticity->execute();
		timerBulk.stop();
		timerBC.resume();
		executeAll(bc);
		timerBC.stop();		
		if (!(i % tout.v()))
		{
			cout << i << endl;
			writer.write();
		}
	}
	timer.stop();
	
	std::cout << "Finished" << endl;	

	cout << "time=" << timer.getTime() << "; clockTime="
		 <<  timer.getClockTime() <<  "; load=" 
		 <<  timer.getProcessorLoad() * 100 << "%" << endl;
	cout << "timeBulk=" << timerBulk.getTime() << 
		    "; timeBC=" << timerBC.getTime() << endl;

	std::cout << "Output...";
	std::cout << "Finished" << endl;	
	std::cout << "Ok" << endl;

	return 0;
}
