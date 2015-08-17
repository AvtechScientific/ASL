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


#include "aclHardware.h"
#include "../aslUtilities.h"
#include "Kernels/aclKernel.h"
#include <fstream>


using namespace std;
using namespace asl;

namespace acl
{
	vector<string> clExtension({"cl_khr_fp64",
								"cl_khr_int64_base_atomics",
								"cl_khr_int64_extended_atomics",
								"CL_KHR_gl_sharing"});

	// order of elements in TYPE and TYPE_SIZE
	// has to correspond to the order of elements of TypeID
	// don't change strings to "cl_..."!!!
	const vector<string> TYPE({"int", "uint", "float", "double", "long"});
	const vector<unsigned char> TYPE_SIZE({sizeof(cl_int), sizeof(cl_uint), sizeof(cl_float), sizeof(cl_double), sizeof(cl_long)});
	const vector<TypeID> TYPE_SELECT({TYPE_INT, TYPE_INT, TYPE_INT, TYPE_LONG, TYPE_LONG});

	Hardware hardware;


	Hardware::Hardware():
		devicesInfo("")
	{
		// Scan all available platforms and their devices
		vector<cl::Platform> platforms;
		vector<cl::Device> devices;
		cl_context_properties cps[3];
		cl::Context context;
		
		cl_int status = 0;	
		status = cl::Platform::get(&platforms);
		errorMessage(status, "acl::Platform::get()");
		
		if (platforms.size() > 0)
		{
			for (unsigned int i = 0; i < platforms.size(); ++i)
			{
				// Probably not needed, since getDevices() seems to call clear() internally
				devices.clear();
				status = platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &devices);
				errorMessage(status, "acl::Platform::getDevices()");
				devicesInfo += "Platform: " + platforms[i].getInfo<CL_PLATFORM_VENDOR>()
					    + "\nNumber of devices: " + numToStr(devices.size()) + "\n";

				cps[0] = CL_CONTEXT_PLATFORM;
				cps[1] = (cl_context_properties)(platforms[i])();
				cps[2] = 0;

				for (unsigned int j = 0; j < devices.size(); ++j)
				{
					// Create an OpenCL context for the current device
					context = cl::Context(vector<cl::Device>(1, devices[j]), cps, NULL, NULL, &status);
					errorMessage(status, "acl::Context::Context()");

					// Create an OpenCL command queue for current context and device
					queues.push_back(CommandQueue(new cl::CommandQueue(context, devices[j], 0, &status)));
					errorMessage(status, "acl::CommandQueue::CommandQueue()");
 					
					devicesInfo += "\t" + devices[j].getInfo<CL_DEVICE_NAME>() + "\n";
				}
			}
		}

		// Set first found queue as default queue
		defaultQueue = queues.front();
	}


	void Hardware::setDefaultQueue(const std::string & platform,
								   const std::string & device)
	{
		defaultQueue = nullptr;

		for (unsigned int i(0); i < queues.size(); ++i)
		{
			if ((platform == getPlatformVendor(queues[i])) && 
			    (device == getDeviceName(queues[i])))
			{
				// Choose requested device on requested platform
				defaultQueue = queues[i];
			}
		}

		// Warn if requested combination of platform and device was not found
		if (defaultQueue.get() == nullptr)
		{
			// Choose first available device
			defaultQueue = queues.front();
			warningMessage("Requested combination of platform(" + platform + ") and device(" + device + ") not found! Using:\n" + getDefaultDeviceInfo());
		}
	}


	std::string Hardware::getDefaultDeviceInfo()
	{
		string defaultDeviceInfo("platform = " + getPlatformVendor(defaultQueue)
								+ "\ndevice = " + getDeviceName(defaultQueue));
		return defaultDeviceInfo;
	}


	string Hardware::getDevicesInfo()
	{
		return devicesInfo;
	}
	

	string getPlatformVendor(const CommandQueue & queue)
	{
		cl_context_properties cps = getContext(queue).getInfo<CL_CONTEXT_PROPERTIES>()[1];
		return  (cl::Platform((cl_platform_id)cps)).getInfo<CL_PLATFORM_VENDOR>();
	}


	string getDeviceName(const CommandQueue & queue)
	{
		return getDevice(queue).getInfo<CL_DEVICE_NAME>();
	}


	string getDeviceVersion(const CommandQueue & queue)
	{
		return getDevice(queue).getInfo<CL_DEVICE_VERSION>();
	}


	cl_device_type getDeviceType(const CommandQueue & queue)
	{
		return getDevice(queue).getInfo<CL_DEVICE_TYPE>();
	}


	cl_uint getNComputeUnits(const CommandQueue & queue)
	{
		return getDevice(queue).getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
	}


	cl::Device getDevice(const CommandQueue & queue)
	{
		return queue->getInfo<CL_QUEUE_DEVICE>();
	}


	cl::Context getContext(const CommandQueue & queue)
	{
		return queue->getInfo<CL_QUEUE_CONTEXT>();
	}


	unsigned int getAlignment(const CommandQueue & queue)
	{
		return getDevice(queue).getInfo<CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE>();
	}


	size_t getMaxItemSize(const CommandQueue & queue)
	{
		return getDevice(queue).getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0];
	}


	cl_uint getVectorWidth(const CommandQueue & queue, const TypeID typeID)
	{
		cl_uint width;

		switch (typeID)
		{
			case TYPE_INT:
				width = getDevice(queue).getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_INT>();
				break;
			case TYPE_UINT:
				// same width as in TYPE_INT
				width = getDevice(queue).getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_INT>();
				break;
			case TYPE_FLOAT:
				width = getDevice(queue).getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT>();
				break;
			case TYPE_LONG:
				width = getDevice(queue).getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG>();
				break;
			default:
				width = getDevice(queue).getInfo<CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE>();
				break;
		}

		return width;
	}
	

	bool extensionAvailable(const CommandQueue & queue, const Extension extension)
	{
		string availableExtensions = getDevice(queue).getInfo<CL_DEVICE_EXTENSIONS>();
		return (availableExtensions.find(clExtension[extension]) != string::npos);	
	}


	cl_device_local_mem_type getLocalMemoryType(const CommandQueue & queue)
	{
		return getDevice(queue).getInfo<CL_DEVICE_LOCAL_MEM_TYPE>();
	}


	cl_ulong getLocalMemorySize(const CommandQueue & queue)
	{
		return getDevice(queue).getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
	}


	cl_ulong getKernelLocalMemSize(const Kernel & kernel)
	{
		return kernel.getKernel().getWorkGroupInfo<CL_KERNEL_LOCAL_MEM_SIZE>(getDevice(kernel.getQueue()));
	}


	cl_ulong getKernelPrivateMemSize(const Kernel & kernel)
	{
		return kernel.getKernel().getWorkGroupInfo<CL_KERNEL_PRIVATE_MEM_SIZE >(getDevice(kernel.getQueue()));
//		return kernel.getKernel().getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE >(getDevice(kernel.getQueue())); 

	}


	cl_device_fp_config doublePrecisionSupport(const CommandQueue & queue)
	{
		return getDevice(queue).getInfo<CL_DEVICE_DOUBLE_FP_CONFIG>();
	}


} // namespace acl
