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


#ifndef ASLNUMMETHODSMERGER_H
#define ASLNUMMETHODSMERGER_H

#include <memory>
#include <vector>

namespace acl
{
	class KernelMerger;
	typedef std::shared_ptr<KernelMerger> SPKernelMerger;
}

namespace asl
{
	class SingleKernelNM;
	typedef std::shared_ptr<SingleKernelNM> SPSingleKernelNM;
	/// Mergers several methods into single kernel
	/// \ingroup Numerics
	class NumMethodsMerger
	{
		public:
			acl::SPKernelMerger km;
		private:
			std::vector<SPSingleKernelNM> nmList;
		public:
			NumMethodsMerger();
			
			/// Executes the numerical procedures
			void execute();
			/// Builds the necesery internal data and kernels
			/** \param buildKernels defines whether to call init0 or not*/
			void init(bool buildKernels);
			
			inline void addNM(SPSingleKernelNM nm);
			void addNM(const std::vector<SPSingleKernelNM> & nm);
	};

//---------------------------- Implementation ---------------------------

	inline void NumMethodsMerger::addNM(SPSingleKernelNM nm)
	{
		nmList.push_back(nm);
	}
	
}	//asl

#endif //ASLNUMMETHODSMERGER_H
