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
	asl-hardware utility.
	Displays hardware information.
 */

#include "acl/aclHardware.h"
#include "aslUtilities.h"

using namespace acl;
using namespace std;
using namespace asl;

string typeToString(unsigned int t)
{
	string s;
	switch (t)
	{
		case CL_DEVICE_TYPE_CPU : s="CPU"; break;
		case CL_DEVICE_TYPE_GPU : s="GPU"; break;
		case CL_DEVICE_TYPE_ACCELERATOR : s="ACCELERATOR"; break;
		case CL_DEVICE_TYPE_DEFAULT : s="DEFAULT"; break;
//		case CL_DEVICE_TYPE_CUSTOM : s="CUSTOM"; break; //in opencl 1.1 is undefined
		default: s="type is unknown";
	}
	return s;
}

void printHardwareInfo(const CommandQueue & queue)
{

	cout << "\t\ttype: " << typeToString(getDeviceType(queue)) << endl;
	cout << "\t\tnumber of compute units: " << getNComputeUnits(queue) << endl;
	cout << "\t\talignment: " << getAlignment(queue) << endl;
	cout << "\t\tlocal memory type: "
		 << (getLocalMemoryType(queue) == CL_LOCAL ? "CL_LOCAL" : "CL_GLOBAL") << endl;
	cout << "\t\tlocal memory size: " << getLocalMemorySize(queue) << endl;
	cout << "\t\tmax item size: " << getMaxItemSize(queue) << endl;
	cout << "\t\tvector width float: " << getVectorWidth(queue, TYPE_FLOAT) << endl;
	cout << "\t\tvector width double: " << getVectorWidth(queue, TYPE_DOUBLE) << endl;
	cout << "\t\textension CL_KHR_FP64: "
		 << extensionAvailable(queue, CL_KHR_FP64) << endl;
	cout << "\t\textension CL_KHR_INT64_EXTENDED_ATOMICS: "
		 << extensionAvailable(queue, CL_KHR_INT64_EXTENDED_ATOMICS) << endl;
	cout << "\t\tsupported OpenCL version: " << getDeviceVersion(queue) << endl;
}


int main()
{

	// Have a look at the available platforms and their devices
	vector<cl::Platform> platforms;
	vector<cl::Device> devices;
	cl_context_properties cps[3];
	cl::Context context;
	CommandQueue queue;
		
	cl_int status = 0;	
	status = cl::Platform::get(&platforms);
	errorMessage(status, "Platform::get()");
		
	if (platforms.size() > 0)
	{
		for (unsigned int i = 0; i < platforms.size(); ++i)
		{
			status = platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &devices);
			errorMessage(status, "Platform::getDevices()");
			cout << "Platform: " << platforms[i].getInfo<CL_PLATFORM_VENDOR>()
				 << "\nNumber of devices: " << devices.size() << endl;

			cps[0] = CL_CONTEXT_PLATFORM;
			cps[1] = (cl_context_properties)(platforms[i])();
			cps[2] = 0;

			for (unsigned int j = 0; j < devices.size(); ++j)
			{
				// Create an OpenCL context for the current device
				context = cl::Context(vector<cl::Device>(1, devices[j]), cps, NULL, NULL, &status);
				errorMessage(status, "Context::Context()");

				// Create an OpenCL command queue for current context and device
				queue = CommandQueue(new cl::CommandQueue(context, devices[j], 0, &status));
				errorMessage(status, "CommandQueue::CommandQueue()");

				cout << "\t" << devices[j].getInfo<CL_DEVICE_NAME>() << endl;
				printHardwareInfo(queue);
				cout << endl;
			}
			cout << endl;
		}
	}

	return 0;
}