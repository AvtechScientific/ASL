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
	\example levelSetBasic.cc
*/

#include <aslDataInc.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslInterfaceTrackingAlg1.h>
#include <utilities/aslTimer.h>
#include <utilities/aslParametersManager.h>
#include <math/aslTemplates.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <acl/aclUtilities.h>
#include <aslGeomInc.h>

typedef float FlT;
//typedef double FlT;
typedef asl::UValue<FlT> Param;
acl::TypeID type(acl::typeToTypeID<FlT>()); 


int main(int argc, char* argv[])
{
	asl::ApplicationParametersManager appParamsManager("levelSetBasic", "1.0");

	asl::Parameter<asl::AVec<int>> size(asl::makeAVec<int>(100, 100, 100), "size", "size");
	asl::Parameter<FlT> dx(1.0, "dx", "dx");
	asl::Parameter<FlT> dt(1.0, "dt", "dt");
	asl::Parameter<asl::AVec<FlT>> v(asl::makeAVec<FlT>(0.0, 0.0, 0.0), "v", "v");
	asl::Parameter<FlT> radius(10.5, "radius", "initial radius");
	

	asl::Parameter<cl_uint> nIterations(100, "nIterations", "Number of iterations");
	asl::Parameter<cl_uint> nItOut(10, "nItOut", "Number of iterations for output");

	appParamsManager.load(argc, argv);
	
	std::cout << "Data initialization... ";

	asl::Block block(size.v(), dx.v());
	auto levelSet(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	auto sphere(generateDFSphere(radius.v(), asl::AVec<>(asl::AVec<FlT>(size.v())*FlT(.5))));
	asl::initData(levelSet, -normalize(sphere,dx.v()));
	
	auto velocity(asl::generateDataContainerConst_SP(block, v.v(), 1u));

	
	asl::WriterVTKXML writer(appParamsManager.getDir() + "levelSetBasic");
	writer.addScalars("levelSet", *levelSet);
	
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization... " << flush;

	auto lsNum(std::make_shared<asl::InterfaceTrackingAlg1>(levelSet,velocity));
	
	lsNum->init();

	std::cout << "Finished" << endl;
	std::cout << "Computing...";
	asl::Timer timer;

	writer.write();
	
	timer.start();
	for (unsigned int i(0); i < nIterations.v(); ++i)
	{
		lsNum->execute();
		if(!(i % nItOut.v()))
			writer.write();		
	}
	timer.stop();
	
	cout << "Finished" << endl;	

	cout << "Computation statistic:" << endl;
	cout << "Real Time = " << timer.realTime() << "; Processor Time = "
		 << timer.processorTime() << "; Processor Load = "
		 << timer.processorLoad() * 100 << "%" << endl;

	return 0;
}
