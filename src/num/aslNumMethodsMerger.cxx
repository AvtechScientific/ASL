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


#include "aslNumMethodsMerger.h"
#include "aslSingleKernelNM.h"
#include <acl/Kernels/aclKernelMerger.h>

using namespace std;

namespace asl
{

	NumMethodsMerger::NumMethodsMerger():
		km(new acl::KernelMerger)
	{
	}

	void NumMethodsMerger::init(bool buildKernels)
	{
		for(unsigned int i(0); i<nmList.size();++i){
			if(buildKernels)
				nmList[i]->init0();
			km->addKernel(nmList[i]->kernel);
		}
		km->setup();
	}

	void NumMethodsMerger::execute()
	{
		km->compute();
		for(unsigned int i(0); i<nmList.size();++i){
			nmList[i]->postProcessing();
		}
	}

	void NumMethodsMerger::addNM(const vector<SPSingleKernelNM> &nm)
	{
		unsigned int s(nmList.size());
		nmList.resize(s+nm.size());
		for(unsigned int i(0);i<nm.size();++i)
			nmList[i+s]=nm[i];
	}
		

} // asl

