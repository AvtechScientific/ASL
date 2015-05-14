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


#include "aclKernel.h"
#include "../aclHardware.h"
#include <acl/aclElementBase.h>
#include <iostream>
#include "../../aslUtilities.h"
#include <algorithm>

using namespace std;
using namespace asl;


namespace acl
{

	unsigned int Kernel::kernelNum = 0;
	
	Kernel::Kernel(const KernelConfiguration kernelConfig_):
		groupsNumber(0),
		kernelConfig(kernelConfig_),
		kernelSource("")
	{
		id = kernelNum;
		++kernelNum;
	}


	cl_uint Kernel::detectVectorWidth()
	{
		vector<cl_uint> widths;

		// get all possible vector widths for this queue
		for (unsigned int i = 0; i < TYPE.size(); ++i)
			widths.push_back(getVectorWidth(queue, (TypeID)i));

		cl_uint absoluteMin = *min_element(widths.begin(), widths.end());
		cl_uint min = *max_element(widths.begin(), widths.end());

		unsigned int i = 0;
		while ((min > absoluteMin) && (i < arguments.size()))
		{
			if (min > widths[arguments[i]->getTypeID()])
				min =  widths[arguments[i]->getTypeID()];
			++i;
		}

		i = 0;
		while ((min > absoluteMin) && (i < localDeclarations.size()))
		{
			if (min > widths[localDeclarations[i]->getTypeID()])
				min =  widths[localDeclarations[i]->getTypeID()];
			++i;
		}

		// Workaround for the case when double vector width is 0.
		// The code still works, however with acceleration
		// only for types other than double.
		min = min == 0 ? 1 : min;
		return min;
	}

			
	// Enables double extension only if there are doubles in the source,
	// enables right extension (workaround for the incomplete AMD cl_amd_fp64)
	void enableDoubleExtension(string & kernelSource, const CommandQueue & queue)
	{
		if (kernelSource.find(TYPE[TYPE_DOUBLE]) != string::npos)
		{
			string availableExtensions = getDevice(queue).getInfo<CL_DEVICE_EXTENSIONS>();

			// workaround for the incomplete AMD cl_amd_fp64.
			// once AMD implements this extension completely
			// only the else-block should stay
			if (availableExtensions.find("cl_amd_fp64") != string::npos)
			{
				kernelSource = "#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n" + kernelSource;
			}
			else
			{
				kernelSource = "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n" + kernelSource;
			}
		}
		else
		{
			string availableExtensions = getDevice(queue).getInfo<CL_DEVICE_EXTENSIONS>();

			// workaround for the incomplete AMD cl_amd_fp64.
			// once AMD implements this extension completely
			// only the else-block should stay
			if (availableExtensions.find("cl_amd_fp64") != string::npos)
			{
				kernelSource = "#pragma OPENCL EXTENSION cl_amd_fp64 : disable\n" + kernelSource;
			}
			else
			{
				kernelSource = "#pragma OPENCL EXTENSION cl_khr_fp64 : disable\n" + kernelSource;
			}
			
		}
	}


	void Kernel::generateExtensions()
	{
		kernelSource = "\n" + kernelSource;
		for (unsigned int i = 0; i < kernelConfig.extensions.size(); i++)
		{
			kernelSource = "#pragma OPENCL EXTENSION " + kernelConfig.extensions[i]
							+ " : enable\n" + kernelSource;
		}
		enableDoubleExtension(kernelSource, queue);
	}

		
	void Kernel::generateArguments()
	{
		for (unsigned int i = 0; i < arguments.size(); i++)
		{
			kernelSource += arguments[i]->getTypeSignature(kernelConfig) + ",\n                       ";
			for (unsigned int j = 0; j < numToStr(id).size(); ++j)
			{
				kernelSource += " ";
			}
		}
		// drop all after the ',' of last parameter
		kernelSource.erase(kernelSource.size() - 25 - numToStr(id).size());
		kernelSource += ")\n{\n\t";
	}


	void Kernel::generateIndex()
	{
		if (kernelConfig.local)
		{
			kernelSource += "uint " + INDEX + " = get_local_id(0);\n\t";
			kernelSource += "uint groupID = get_group_id(0);\n";
		}
		else
		{
			if ((kernelConfig.vectorWidth > 1) && (kernelConfig.unaligned))
			{
				kernelSource += "uint " + INDEX + " = "
							 + numToStr(kernelConfig.vectorWidth)
							 + " * get_global_id(0);\n";
			}
			else
			{
				kernelSource += "uint " + INDEX + " = get_global_id(0);\n";
			}
		}
	}


	void Kernel::generateLocalDeclarations()
	{
		for (unsigned int i = 0; i < localDeclarations.size(); i++)
		{
			kernelSource += "\t" + localDeclarations[i]->getLocalDeclaration(kernelConfig) + ";\n";
		}
	}


	void Kernel::generateExpressions()
	{
		for (unsigned int i = 0; i < expression.size(); i++)
		{
			kernelSource += "\t" + expression[i]->str(kernelConfig)  + ";\n";
		}

		kernelSource += "}";
	}


	void Kernel::generateKernelSource()
	{
		kernelSource = "__kernel void compute_" + numToStr(id) + "(";
		
		filterDeclarations();

		generateArguments();

		generateIndex();

		generateLocalDeclarations();
		
		generateExpressions();

		// enable extensions in the last step after
		// the kernel source is available for analysis
		generateExtensions();
		
		regenerateKernelSource = false;
	}


	void Kernel::buildKernel()
	{
		cl_int status = 0;
		
		cl::Program::Sources source(1, make_pair(kernelSource.data(), kernelSource.size()));

		cl::Program program = cl::Program(getContext(queue), source, &status);
		errorMessage(status, "Program::Program()");
		status = program.build(std::vector<cl::Device>(1, getDevice(queue)));

		if (status == CL_BUILD_PROGRAM_FAILURE)
		{
			string str = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(getDevice(queue));

			cout << " \n\t\t\tBUILD LOG\n";
			cout << " ************************************************\n";
			cout << str << endl;
			cout << " ************************************************\n";
			cout << " \n\t\t\tKERNEL SOURCE CODE\n";
			cout << " ------------------------------------------------\n";
			cout << kernelSource << endl;
			cout << " ------------------------------------------------\n";
		}

		errorMessage(status, "Program::build()");

		string kernelName;
		kernelName = "compute_" + numToStr(id);
		kernel = cl::Kernel(program, kernelName.c_str(), &status);
		errorMessage(status, "Kernel() - kernelBase_" + numToStr(id));
//		cout<<"  "<<kernelName<<": privMem="<< getKernelPrivateMemSize(*this)<<endl;
		
	}


	void Kernel::setKernelArguments()
	{
		for (unsigned int i = 0; i < arguments.size(); ++i)
		{
			arguments[i]->setAsArgument(kernel, i);
		}
	}


	void Kernel::updateKernelConfiguration()
	{
		// if it is a simd kernel - detect proper vector width
		if (kernelConfig.vectorWidth > 1) 
			kernelConfig.vectorWidth = detectVectorWidth();
	}


	void Kernel::setup()
	{
		if (size == 0)
			errorMessage("Kernel::setup() - kernel's size is 0; add proper expression");

		if (kernelConfig.local && (size > getMaxItemSize(queue)))
			errorMessage("Kernel::setup() - requested group size is larger than supported by device");

		if (kernelConfig.local && (groupsNumber == 0))
			errorMessage("Kernel::setup() - groups number was not set");

		updateKernelConfiguration();
		generateKernelSource();
		buildKernel();
	}


	inline unsigned int paddingElementsSIMD(unsigned int size, unsigned int vectorWidth);
	
	void Kernel::compute()
	{
		if (regenerateKernelSource)
			setup();
	
		cl_int status = 0;
		cl::Event event;

		setKernelArguments();

		size_t paddedSize = (size + paddingElementsSIMD(size, kernelConfig.vectorWidth)) / kernelConfig.vectorWidth;
		// if 4th argument is cl::NullRange - OpenCL implementation
		// will determine how to break the global work-items into
		// appropriate work-group instances
		status = queue->enqueueNDRangeKernel(kernel,
	 	                                     cl::NullRange,
	    	                                 kernelConfig.local ? cl::NDRange(groupsNumber * paddedSize) : cl::NDRange(paddedSize),
	    	                                 kernelConfig.local ? cl::NDRange(paddedSize) : cl::NullRange,
	    	                                 NULL,
	    	                                 &event);
		errorMessage(status, "CommandQueue::enqueueNDRangeKernel() - kernel");

		status = event.wait();
		errorMessage(status, "Event::wait() - event");
	}


	void Kernel::setGroupsNumber(unsigned int n)
	{
		groupsNumber = n;
	}

	unsigned int Kernel::getGroupsNumber()
	{
		return groupsNumber;
	}		

	string Kernel::getKernelSource()
	{
		return kernelSource;
	}


	unsigned int Kernel::getKernelID()
	{
		return id;
	}


	const cl::Kernel & Kernel::getKernel() const
	{
		return kernel;
	}


	void Kernel::clear()
	{
		kernelSource = "__kernel void compute_" + numToStr(id) + "(";
		regenerateKernelSource = true;
		expression.clear();
		localDeclarations.clear();
		arguments.clear();
		size = 0;
	}


	inline unsigned int paddingElementsSIMD(unsigned int size, unsigned int vectorWidth)
	{
		// second modulo is added in order to make padding = 0 in the case
		// that size is divisible by vectorWidth
		return (vectorWidth - (size % vectorWidth)) % vectorWidth;
	}
		
} // namespace acl
