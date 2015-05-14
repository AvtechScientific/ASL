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
 	\example testPerformance.cc

	\todo repare
 
	Results:
 	processor: Intel(R) Core(TM)2 Duo CPU     T8100  @ 2.10GHz
	\code
	Test of Sum Operator...
	GPU: 9
	CPU: 12
	\endcode
*/

#include "acl/acl.h"
#include "acl/aclHardware.h"
#include "acl/DataTypes/aclIndex.h"
#include "acl/DataTypes/aclConstant.h"
#include "acl/DataTypes/aclVariable.h"
#include "acl/DataTypes/aclPrivateVariable.h"
#include "acl/DataTypes/aclArray.h"
#include "acl/DataTypes/aclSubvector.h"
#include "acl/Kernels/aclKernel.h" 
#include "aslUtilities.h"
#include "utilities/aslTimer.h"
#include <math.h>
#include <fstream>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>

#define ARRAY_SIZE 10000000
//#define ITERATIONS_NUM 20
#define TIME_INTERVAL 5 // in seconds

using namespace acl;
using namespace std;


// Operators: +
template <typename T> inline T testSum(T x1, T x2)
{ 
	using namespace elementOperators;
	return x1 + x2;
}


// Operators: + (non sequential)
template <typename T> inline T testSumNonSequential(T x1, T x2)
{
	Element index(new Index(ARRAY_SIZE));
	Element c(new Constant<cl_uint>(17));
	Element cArraySize(new Constant<cl_uint>(ARRAY_SIZE));
	using namespace elementOperators;
	return excerpt(x1 + x2, (index * c) % cArraySize);
}


// Operators: +, -, *, /
template <typename T> inline T testBasicOperators(T x1, T x2)
{
	using namespace elementOperators;
	return x1 * (x2 + x2) * x1 + (x2 + x1 * x1) * x2 * x1 * x2 * x1 - x1 / x2 + 
		    x1 * (x1 + x2) * x1 + (x2 - x1 * x2) * x2 - x1 * x2 * x1 - x2 / x1 + 
		    x1 * (x2 + x1) * x1 - (x2 + x1 * x1) * x2 * x1 * x2 * x1 + x1 / x2 -
		    x1 * (x1 + x2) * x2 + (x2 - x1 * x2) * x2 - x1 * x2 * x1 - x2 / x1 + 
		    x1 * (x2 + x1) * x1 - (x2 + x1 * x1) * x2 * x1 * x2 * x1 + x1 / x2;
}


// Operators: +, -, *, /, sin, cos, sqrt
template <typename T> inline T testSpecialOperators(T x1, T x2)
{
	using namespace elementOperators;
	return cos(x1*sin(x2+x2)*x1+sqrt(x2+x1*x1)*x2*sin(x1)*x1*sqrt(x2)*x2*x1*sqrt(x1*x2)-x1/x2);
}


/*
void print(Glib::KeyFile * keyFile,
           string deviceName,
           string key,
           string value,
           bool writeToKeyFile)
{
	cout << key + ": " + value << endl;
	
	// Save the output for this device to the performance report file if needed
	if (writeToKeyFile)
	{
		keyFile->set_string((deviceName),
		                    key,
		                    value);
	}
}
*/

/*
template <typename T> inline void testKernelPerformance(const KernelConfiguration kernelConfig,
                                                        const CommandQueue & queue,
                                                        Glib::KeyFile * keyFile,
                                                        bool writeToKeyFile)
{
	Element arr1(new Array<T>(ARRAY_SIZE, queue));
	Element arr2(new Array<T>(ARRAY_SIZE, queue));
	Element arrResult(new Array<T>(ARRAY_SIZE, queue));
	Element c(new Constant<T>(5));

	string typeStr = typeToStr<T>();
	
	string kernelConfigStr;
	if (kernelConfig == KERNEL_BASIC)
		kernelConfigStr = "KERNEL_BASIC";
	else
		if (kernelConfig == KERNEL_SIMD)
			kernelConfigStr = "KERNEL_SIMD";
		else
			if (kernelConfig == KERNEL_SIMDUA)
				kernelConfigStr = "KERNEL_SIMDUA";

	
	Kernel k(kernelConfig);
	
	// Initialization of the arrays
	{
		using namespace elementOperators;
		k.addExpression(operatorAssignmentSafe(arr1, c));
		k.addExpression(operatorAssignmentSafe(arr2, c));
	}
	k.compute();

	// Test of sum
	k.clear();
	{ 
		using namespace elementOperators;
		k.addExpression(operatorAssignmentSafe(arrResult, testSum(arr1, arr2)));
	}
	k.setup();
	
	asl::Timer timer;
	unsigned int iterationsNum;

	timer.start();
	timer.stop();
	iterationsNum = 0;
	do
	{
		timer.resume();
		k.compute();
		++iterationsNum;
		timer.stop();
	}
	while (timer.getTime() < TIME_INTERVAL);

	print(keyFile,
	      getDeviceName(queue),
	      kernelConfigStr + "---" + typeStr + "---" + "Sum",
	      numToStr((float) iterationsNum / timer.getTime()),
	      writeToKeyFile);

	// Test of non sequential sum (only for non SIMD kernels)
	if (kernelConfig.vectorWidth == 1)
	{
		k.clear();
		{ 
			using namespace elementOperators;
			k.addExpression(operatorAssignmentSafe(arrResult,
			                                       testSumNonSequential(arr1, arr2)));
		}
		k.setup();
		
		timer.start();
		timer.stop();
		iterationsNum = 0;
		do
		{
			timer.resume();
			k.compute();
			++iterationsNum;
			timer.stop();
		}
		while (timer.getTime() < TIME_INTERVAL);

		print(keyFile,
		      getDeviceName(queue),
		      kernelConfigStr + "---" + typeStr + "---" + "NonSequentialSum",
		      numToStr((float) iterationsNum / timer.getTime()),
		      writeToKeyFile);

	}

	// Test of basic operators
	k.clear();
	{ 
		using namespace elementOperators;
		k.addExpression(operatorAssignmentSafe(arrResult,
		                                       testBasicOperators(arr1, arr2)));
	}
	k.setup();
	
	timer.start();
	timer.stop();
	iterationsNum = 0;
	do
	{
		timer.resume();
		k.compute();
		++iterationsNum;
		timer.stop();
	}
	while (timer.getTime() < TIME_INTERVAL);

	print(keyFile,
	      getDeviceName(queue),
	      kernelConfigStr + "---" + typeStr + "---" + "BasicOperators",
	      numToStr((float) iterationsNum / timer.getTime()),
	      writeToKeyFile);


	// Test of special operators
	k.clear();
	if (arr1->getTypeID() != TYPE_INT)
	{
		{ 
			using namespace elementOperators;
			k.addExpression(operatorAssignmentSafe(arrResult,
			                                       testSpecialOperators(arr1, arr2)));
		}
		k.setup();
	
		timer.start();
		timer.stop();
		iterationsNum = 0;
		do
		{
			timer.resume();
			k.compute();
			++iterationsNum;
			timer.stop();
		}
		while (timer.getTime() < TIME_INTERVAL);

		print(keyFile,
		      getDeviceName(queue),
		      kernelConfigStr + "---" + typeStr + "---" + "SpecialOperators",
		      numToStr((float) iterationsNum / timer.getTime()),
		      writeToKeyFile);
	}
}
*/

int main()
{
	/*
	const string FILE_NAME("./PerformanceReport.txt");
	bool writeToKeyFile = false;
	bool dumpKeyFile = false;
	
	Glib::KeyFile *keyFile = new Glib::KeyFile;

	ifstream fileCheck(FILE_NAME);
	if (fileCheck.good())
	{
 		keyFile->load_from_file(FILE_NAME);
	}
	else
	{
		warningMessage("Failed to open " + FILE_NAME + " . Creating new file");
		keyFile->set_comment("\n Performance test report\n\n");
	}

	
	for (unsigned int i = 0; i < hardware.queues.size(); ++i)
	{
		if (!keyFile->has_group(getDeviceName(hardware.queues[i])))
		{
			writeToKeyFile = true;
			dumpKeyFile = true;
		}
		
		cout << "\nDevice: " + getDeviceName(hardware.queues[i]) + "\n" << endl;
		
		testKernelPerformance<cl_int>(KERNEL_BASIC,
		                              hardware.queues[i],
		                              keyFile,
		                              writeToKeyFile);
		testKernelPerformance<cl_float>(KERNEL_BASIC,
		                                hardware.queues[i],
		                                keyFile,
		                                writeToKeyFile);
		if (doublePrecisionSupport(hardware.queues[i]))
			testKernelPerformance<cl_double>(KERNEL_BASIC,
			                                 hardware.queues[i],
			                                 keyFile,
			                                 writeToKeyFile);


		testKernelPerformance<cl_int>(KERNEL_SIMD,
		                              hardware.queues[i],
		                              keyFile,
		                              writeToKeyFile);
		testKernelPerformance<cl_float>(KERNEL_SIMD,
		                                hardware.queues[i],
		                                keyFile,
		                                writeToKeyFile);
		if (doublePrecisionSupport(hardware.queues[i]))
			testKernelPerformance<cl_double>(KERNEL_SIMD,
			                                 hardware.queues[i],
			                                 keyFile,
			                                 writeToKeyFile);


		testKernelPerformance<cl_int>(KERNEL_SIMDUA,
		                              hardware.queues[i],
		                              keyFile,
		                              writeToKeyFile);
		testKernelPerformance<cl_float>(KERNEL_SIMDUA,
		                                hardware.queues[i],
		                                keyFile,
		                                writeToKeyFile);
		if (doublePrecisionSupport(hardware.queues[i]))
			testKernelPerformance<cl_double>(KERNEL_SIMDUA,
			                                 hardware.queues[i],
			                                 keyFile,
			                                 writeToKeyFile);

		writeToKeyFile = false;
	}

	if (dumpKeyFile)
	{
		ofstream file;
		file.open(FILE_NAME, ios::app);
		if (!file.good())
		{
			errorMessage("Opening file " + FILE_NAME + " failed");
		}
		else
		{
			file << keyFile->to_data();
			if (!file.good())
				errorMessage("Writing to file " + FILE_NAME + " failed");
		}
		file.close();
	}

	delete keyFile;
	*/
	return 0;
}
