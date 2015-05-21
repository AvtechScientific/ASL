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


#include "aslDFOptimizer.h"
#include <aslGenerators.h>
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include "acl/aclHardware.h"
#include <math/aslTemplateVE.h>
#include <acl/aclMath/aclVectorOfElements.h>

#include <data/aslDataWithGhostNodes.h>
#include <acl/Kernels/aclKernel.h>
#include <math/aslTemplates.h>
#include <math/aslDistanceFunctionAlg.h>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>

namespace asl
{
	DFOptimizer::DFOptimizer():
//		SingleKernelNM(acl::KERNEL_SIMDUA),
		SingleKernelNM(acl::KERNEL_BASIC),
		vectorTemplate(NULL)
	{
	}


	DFOptimizer::DFOptimizer(Data inD, const VectorTemplate* vT):
//		SingleKernelNM(acl::KERNEL_SIMDUA),
		SingleKernelNM(acl::KERNEL_BASIC),
		inData(inD),
		vectorTemplate(vT)
	{
	}
		
	void DFOptimizer::init0()
	{
		unsigned int nDir(vectorTemplate->vectors.size());
		acl::TypeID type(getElementType(inData->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		auto dataX(generateDCFullSafe(inData, -.9));
		
		acl::VectorOfElements vnew(acl::generateVEPrivateVariable(1, type));
		acl::VectorOfElements isAnyComp(acl::generateVEPrivateVariable(1, typeI));
		acl::VectorOfElements isAnyGhost(acl::generateVEPrivateVariable(1, typeI));
		
		TemplateVE dT;
		dT.init(*dataX, *vectorTemplate);
		(*kernel) << dT.initValues;
		(*kernel) << (isAnyComp = acl::generateVEConstant(0));
		(*kernel) << (isAnyGhost = acl::generateVEConstant(0));
		for(unsigned int i(1); i < nDir; ++i)
		{
			(*kernel) << (isAnyComp = isAnyComp || isComputationNode(dT, i));
			(*kernel) << (isAnyGhost = isAnyGhost || isGhostNode(dT, i));
		}
		(*kernel) << (vnew = dT.getValue(0));
		(*kernel) << (vnew = select(vnew, acl::generateVEConstant(-0.9998), 
		                            (dT.getValue(0)<=-0.9999) && isAnyComp, type));
		(*kernel) << (vnew = select(vnew, acl::generateVEConstant(0.9998), 
		                            dT.getValue(0)>=0.9999 && isAnyGhost, type));
		(*kernel) << (vnew = select(vnew, acl::generateVEConstant(-1), 
		                            (dT.getValue(0)<0 && dT.getValue(0)>-0.9999) && !isAnyComp,type));
		(*kernel) << (vnew = select(vnew, acl::generateVEConstant(1), 
		                            (dT.getValue(0)>0 && dT.getValue(0)<0.9999) && !isAnyGhost,type));
		
		(*kernel) << (assignmentSafe(inData->getEContainer(), vnew));		
	}
			
	void optimizeMap(SPDataWithGhostNodesACLData c, const VectorTemplate* vt)
	{
		auto nm(make_shared<DFOptimizer> (c, vt));
		nm->init();
		nm->execute();
	}
		
		
} //asl
