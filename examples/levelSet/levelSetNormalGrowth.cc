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
	\example levelSetNormalGrowth.cc
 */

#include <aslDataInc.h>
#include <writers/aslVTKFormatWriters.h>
#include <num/aslLSNormalGrowth.h>
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
	asl::ApplicationParametersManager appParamsManager("levelSetNormalGrowth",
	                                                   "1.0");

	asl::Parameter<asl::AVec<int>> size("size", "size");
	asl::Parameter<FlT> dx("dx", "dx");
	asl::Parameter<FlT> dt("dt", "dt");
	asl::Parameter<FlT> superS("superS", "super satuation");
	asl::Parameter<FlT> radius("radius", "initial radius");
	

	asl::Parameter<cl_uint> nIterations("nIterations", "Number of iterations");
	asl::Parameter<cl_uint> nItOut("nItOut", "Number of iterations for output");

	appParamsManager.load(argc, argv);
	
	std::cout << "Data initialization...";

	asl::Block block(size.v(), dx.v());
	auto levelSet(asl::generateDataContainerACL_SP<FlT>(block, 1, 1u));

	asl::AVec<> center(asl::AVec<FlT>(size.v())*FlT(.5));
		
	auto sphere1(generateDFSphere(radius.v(), center*.8));
	auto sphere2(generateDFSphere(radius.v(), center*1.2));
	asl::initData(levelSet, normalize(sphere1 | sphere2, dx.v()));
	
	auto superSaturation(asl::generateDataContainerConst_SP(block, superS.v(), 1u));

	
	asl::WriterVTKXML writer(appParamsManager.getDir() + "levelSetNormalGrowth");
	writer.addScalars("levelSet", *levelSet);
	
	std::cout << "Finished" << endl;
	
	std::cout << "Numerics initialization..." << flush;

	auto lsNum(std::make_shared<asl::LSNormalGrowth>(levelSet, superSaturation));
	
	lsNum->init();

	std::cout << "Finished" << endl;
	std::cout << "Computing...";
	asl::Timer timer;

	writer.write();
	
	timer.start();
	for (unsigned int i(0); i < nIterations.v(); ++i)
	{
		lsNum->execute();
		if (!(i % nItOut.v()))
			writer.write();		
	}
	timer.stop();
	
	std::cout << "Finished" << endl;	

	cout << "time=" << timer.getTime() << "; clockTime="
		 << timer.getClockTime()	 << "; load=" 
		 << timer.getProcessorLoad() * 100 << "%" << endl;

	std::cout << "Output...";
	std::cout << "Finished" << endl;	
	std::cout << "Ok" << endl;

	return 0;
}
