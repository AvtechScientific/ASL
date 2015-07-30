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


#ifndef ACLHARDWARE_H
#define ACLHARDWARE_H


//#include <CL/cl.hpp>
// Supply "cl.hpp" with ASL, since it is not present in OpenCL 2.0
// Remove the file after switching to OpenCL 2.1
#include "cl.hpp"
#include "aclStdIncludes.h"
#include <memory>
#include "aclTypes.h"

namespace acl
{

	typedef std::shared_ptr<cl::CommandQueue> CommandQueue;


	extern const std::vector<std::string> TYPE;
	extern const std::vector<unsigned char> TYPE_SIZE;
	/// contains trasnlation of types necessery for use in the function select
	/* expression in the condition field should have this type*/ 
	extern const std::vector<TypeID> TYPE_SELECT;

	/// Returns vendor name.
	/// \ingroup HardwareInformation
	std::string getPlatformVendor(const CommandQueue & queue);

	/// Returns device name.
	/// \ingroup HardwareInformation
	std::string getDeviceName(const CommandQueue & queue);

	/** Returns the OpenCL version supported by the device */
	std::string getDeviceVersion(const CommandQueue & queue);

	/// Returns device type.
	/// \ingroup HardwareInformation
	cl_device_type getDeviceType(const CommandQueue & queue);

	/// Returns number of computer units on the device.
	/// \ingroup HardwareInformation
	cl_uint getNComputeUnits(const CommandQueue & queue);
	
	cl::Device getDevice(const CommandQueue & queue);
	cl::Context getContext(const CommandQueue & queue);

	/// Returns the smallest alignment in bytes which can be used for any data type.
	/// \ingroup HardwareInformation
	unsigned int getAlignment(const CommandQueue & queue);
	

	/// Returns type of local memory supported.
	/// CL_LOCAL (implies dedicated local memory storage such as SRAM) or CL_GLOBAL.
	/// \ingroup HardwareInformation
	cl_device_local_mem_type getLocalMemoryType(const CommandQueue & queue);


	/// Size of local memory arena in bytes. The minimum value is 1 KB.
	/// \ingroup HardwareInformation
	cl_ulong getLocalMemorySize(const CommandQueue & queue);


	/// Returns the maximum number of work-items that can be specified
	/// in '0' dimension of the work-group to clEnqueueNDRangeKernel.
	/// \ingroup HardwareInformation
	size_t getMaxItemSize(const CommandQueue & queue);


	/// Returns the native ISA vector width. The vector width is defined as the
	/// number of scalar elements that can be stored in the vector.
	/// If the cl_khr_fp64 extension is not supported,
	/// CL_DEVICE_NATIVE_DOUBLE must return 0.
	/// \ingroup HardwareInformation
	cl_uint getVectorWidth(const CommandQueue & queue, const TypeID typeID);


	/// Checks availability of an OpenCL extension 
	/// \ingroup HardwareInformation
	bool extensionAvailable(const CommandQueue & queue, const Extension extension);


	/// Describes double precision floating-point capability of the OpenCL device.
	/// Returns a non-zero value if double precision FP is supported.
	/// See CL_DEVICE_DOUBLE_FP_CONFIG for more info.
	/// \ingroup HardwareInformation
	cl_device_fp_config doublePrecisionSupport(const CommandQueue & queue);

	class Kernel;
	
	/**
	Returns the amount of local memory in bytes being used by a kernel. This
	includes local memory that may be needed by an implementation to execute
	the kernel, variables declared inside the kernel with the __local address
	qualifier and local memory to be allocated for arguments to the kernel
	declared as pointers with the __local address qualifier and whose size is
	specified with clSetKernelArg.

	\ingroup HardwareInformation
	*/
	cl_ulong getKernelLocalMemSize(const Kernel & kernel);


	/**
	Returns the minimum amount of private memory, in bytes, used by each work-item
	in the kernel. This value may include any private memory needed by
	an implementation to execute the kernel, including that used by the language
	built-ins and variable declared inside the kernel with the __private qualifier.

	\ingroup HardwareInformation
	*/
	cl_ulong getKernelPrivateMemSize(const Kernel & kernel);


	/// Provides access to the underlying hardware
	class Hardware
	{
		public:
			/// OpenCL related initializations are done here.
			/// Context, Device list, Command Queue are set up.
			/// Default computation queue is set to the first found device.
			Hardware();
			/// Sets default computation queue
			/// identified by \p platform and \p device.
			/// Warns if requested combination is not found.
			void setDefaultQueue(const std::string & platform = "",
			                     const std::string & device = "");
			std::vector<CommandQueue> queues;
			CommandQueue defaultQueue;
			std::string getDevicesInfo();
			std::string getDefaultDeviceInfo();
		private:
			std::string devicesInfo;
	};


	// GLOBALS
	extern Hardware hardware;

} // namespace acl
#endif // ACLHARDWARE_H
