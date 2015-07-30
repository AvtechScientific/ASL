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


#ifndef ACLKERNEL_H
#define ACLKERNEL_H

#include "aclExpressionContainer.h"
#include "aclKernelConfiguration.h"
//#include <CL/cl.hpp>
// Supply "cl.hpp" with ASL, since it is not present in OpenCL 2.0
// Remove the file after switching to OpenCL 2.1
#include "acl/cl.hpp"

namespace acl
{

	extern const KernelConfiguration KERNEL_BASIC;
	
	/// OpenCl Kernel generator 
	/**	
		\ingroup KernelGen
		\listedNote The kernel can be run without updating of parameters. 
					This possibility can lead to some optimization.
	 				to realize this possibility the function
		 			\p computeWithoutUpdate can be added

	*/
	class Kernel: public ExpressionContainer
	{
		protected:
			static unsigned int kernelNum;
			unsigned int id;
			unsigned int groupsNumber;
			KernelConfiguration kernelConfig;
			std::string kernelSource;
			cl::Kernel kernel;

			/// detects minimal vector width of all available types of Elements
			cl_uint detectVectorWidth();
			void generateExtensions();
			void generateArguments();
			void generateIndex();
			void generateLocalDeclarations();
			void generateExpressions();
			virtual void generateKernelSource();
			void updateKernelConfiguration();
			void buildKernel();
			void setKernelArguments();
		public:
			explicit Kernel(const KernelConfiguration kernelConfig_ = KERNEL_BASIC);
			/// Prepares kernel for launch.
			/// Should always be called before compute() after all expressions are added.
			/// Generates kernel source, builds kernel and sets its arguments.
			void setup();
			void compute();
			void setGroupsNumber(unsigned int n);
			unsigned int getGroupsNumber();			
			std::string getKernelSource();
			unsigned int getKernelID();
			const cl::Kernel & getKernel() const;
			/// removes all expressions from the kernel
			void clear();
			inline const KernelConfiguration & getConfiguration()const;

			friend class KernelMerger;
	};

	/// \related Kernel
	typedef std::shared_ptr<Kernel> SPKernel;

	/// creates \p n kernels in \p vk with configuration \p kernelConfig_
	/// \related Kernel
	void inline createKernels(std::vector<SPKernel> & vk, 
	                          unsigned int n, 
	                          const KernelConfiguration & kernelConfig_ = KERNEL_BASIC);
		

//---------------------------- Implementation -----------------------------

	void inline createKernels(std::vector<SPKernel> & vk, 
	                          unsigned int n, 
	                          const KernelConfiguration & kernelConfig_)
	{
		vk.resize(n);
		for (unsigned int i(0); i < n; ++i)
		{
			vk[i].reset(new Kernel(kernelConfig_));
		}
	}

	inline const KernelConfiguration & Kernel::getConfiguration()const
	{
		return kernelConfig;
	}
	

} // namespace acl

#endif // ACLKERNEL_H
