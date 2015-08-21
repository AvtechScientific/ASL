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
	\example poroelastic.cc
	Required input file: [brain.vti](http://asl.org.il/input_data/brain.vti)
 */

#include <aslDataInc.h>
#include <acl/aclGenerators.h>
#include <num/aslFDPoroElasticity.h>
#include <num/aslFDElasticityBC.h>
#include <num/aslFDPoroElasticityBC.h>
#include <utilities/aslTimer.h>
#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <aslGeomInc.h>
#include <math/aslDistanceFunction.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <acl/aclUtilities.h>
#include <math/aslIndex2Position.h>
#include <readers/aslVTKFormatReaders.h>
#include <writers/aslVTKFormatWriters.h>


typedef float FlT;
//typedef double FlT;
typedef asl::UValue<FlT> Param;


int main(int argc, char* argv[])
{
	asl::ApplicationParametersManager appParamsManager("poroelastic", "1.0");
	asl::Parameter<asl::AVec<int> > size("size", "size 3D");
	asl::Parameter<cl_float> dx("dx", "dx");
	asl::Parameter<cl_float> dt("dt", "dt");
	asl::Parameter<cl_float> bulkModulus("bulk_modulus", "bulk modulus");
	asl::Parameter<cl_float> shearModulus("shear_modulus", "shear modulus");
	asl::Parameter<cl_float> hydraulicConductivity("hydraulic_conductivity", "hydraulic conductivity");
	asl::Parameter<cl_float> rho("rho", "density");
	asl::Parameter<asl::AVec<FlT> > g("g", "gravity vector");
	asl::Parameter<string> input("input", "path to the brain geometry input file");
	asl::Parameter<unsigned int> tsim("num_iterations", "number of iterations");
	asl::Parameter<unsigned int> tout("num_it_out", "number of iterations between outputs");
	
	appParamsManager.load(argc, argv);
		
	std::cout << "Data initialization... " << flush;

	asl::SPDataWithGhostNodesACLData map0(asl::read(input.v(), 0));
//	asl::Block block(size.v(), dx.v());
	asl::Block block(map0->getInternalBlock());
	dx.v() = block.dx*1e-3;

	dt.v() = 0.0001f;
	Param bulkModulusNum(bulkModulus.v()/rho.v()/dx.v()/dx.v());
	Param shearModulusNum(shearModulus.v()/rho.v()/dx.v()/dx.v());
	asl::AVec<FlT> gNum(g.v()/dx.v());
	Param hydraulicConductivityNum(hydraulicConductivity.v()/dx.v()/dx.v()/dx.v());

	cout << gNum << "; " << bulkModulusNum.v() << "; " << hydraulicConductivityNum.v() << endl;

	
	auto displacement(asl::generateDataContainerACL_SP<FlT>(block, 3, 1u));
	auto pressure(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
	acl::initData(displacement->getEContainer(), acl::generateVEConstantN(3,0));
	acl::initData(pressure->getEContainer(), acl::generateVEConstant(0));
	
//	auto mapDF(asl::generateDFInBlock(block, 0));
//	auto map(asl::generateDataContainer_SP(block, mapDF, 1u, acl::typeToTypeID<FlT>()));
	auto mapX(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));
//	initData(mapX->getEContainer(), map->getEContainer());
	initData(mapX->getEContainer(), map0->getEContainer()*2.-1., acl::KERNEL_BASIC);
	
	asl::WriterVTKXML writer(appParamsManager.getDir() + "poroelastic");
	writer.addVector("displacement", *displacement);
	writer.addScalars("pressure", *pressure);
	writer.addScalars("map", *mapX);
//	writer.write();

	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... " << flush;

	auto elasticity(make_shared<asl::FDPoroElasticity>(displacement,
	                                                   pressure,
   	                                                   acl::generateVEConstant(bulkModulusNum.v()),
	                                                   acl::generateVEConstant(shearModulusNum.v()),
	                                                   acl::generateVEConstant(hydraulicConductivityNum.v()),
	                                                   &asl::d3q15()));

	asl::Index2PositionACL i2p(displacement->getBlock(),acl::typeToTypeID<FlT>());
	asl::AVec<> center(asl::AVec<>(displacement->getBlock().getSize())*.5 *
	                   displacement->getBlock().dx +
	                   displacement->getBlock().position);
	asl::SPDistanceFunction scf(asl::generateDFPlane(asl::AVec<>(g.v()),center));
	
	auto force((1.-acl::sign(scf->getDistance(i2p.positionWithInit)))/2.*
	            acl::generateVEConstant(gNum));
	auto forceField(asl::generateDataContainer_SP(block, force, 1u));
	
	elasticity->setForce(forceField->getSubContainer());
	elasticity->init();

	vector<asl::SPNumMethod> bcl;
	bcl.push_back(asl::generateBCZeroStress(elasticity, mapX));
	asl::addBCRigidWall(bcl, elasticity, {asl::X0});

//!!!!!!	addSliceZ(*bcRigidWall, block.getSize()[0]/2 + 1); return

	initAll(bcl);

	std::cout << "Finished" << endl;
	std::cout << "Computing..." << flush;
	asl::Timer timer;

	executeAll(bcl);

	writer.write();

	timer.start();
	for (unsigned int i(0); i < tsim.v(); ++i)
	{
		elasticity->execute();
		executeAll(bcl);
		if (!(i % tout.v()))
			writer.write();
	}
	timer.stop();
	
	std::cout << "Finished" << endl;	

	cout << "time: " << timer.getTime() << "; clockTime: "
		 <<  timer.getClockTime() <<  "; load: " 
		 <<  timer.getProcessorLoad() * 100 << "%" << endl;

	std::cout << "Output...";
	std::cout << "Finished" << endl;	
	std::cout << "Ok" << endl;

	return 0;
}