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


#ifndef ACLKERNELMERGER_H
#define ACLKERNELMERGER_H

#include <memory>
#include <vector>

namespace acl
{
	class ElementBase;
	typedef std::shared_ptr<ElementBase> Element;
	class Kernel;
	typedef std::shared_ptr<Kernel> SPKernel;

	/// OpenCl Kernel merger 
	/**	
		\ingroup KernelGen
		The KernelMerger generates a big kernel containing differen 
		kernels. The number of copies is a sum of all kernels and 
		it realizes tree like "if" "else" structure in order to execute 
		different kernels
	*/
	class KernelMerger
	{
		private:
			std::vector<SPKernel> kernels;
			SPKernel kernel;
			unsigned int size;
			std::vector<unsigned int> offsets;

			unsigned int getKernelSize(unsigned int i);
			void checkKernelsCompatibility();
			void computeOffsets();
			Element castSpliter(unsigned int i1, unsigned int i2);
		public:
			KernelMerger();
			void setup();
			void compute();
			std::string getKernelSource();
			/// removes all kernels
			void clear();
			void addKernel(SPKernel k);
			void addKernel(const KernelMerger & km);
			unsigned int getSize();
	};

	/// \related Kernel
	typedef std::shared_ptr<KernelMerger> SPKernelMerger;


	
} // namespace acl

#endif // ACLKERNELMERGER_H
