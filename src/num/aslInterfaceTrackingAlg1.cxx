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


#include "aslInterfaceTrackingAlg1.h"
#include <math/aslTemplates.h>
#include <data/aslDataWithGhostNodes.h>
#include "acl/acl.h"
#include "acl/Kernels/aclKernel.h"
#include "acl/Kernels/aclKernelConfigurationTemplates.h"
#include "acl/Kernels/aclKernelMerger.h"
#include "acl/aclGenerators.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include "acl/aclMath/aclBarycentric.h"
#include <algorithm>
#include <math/aslTemplates.h>
#include <math/aslTemplateVE.h>
#include <math/aslTemplatesExtras.h>
#include <math/aslDistanceFunction.h>

namespace asl
{

	InterfaceTrackingAlg1::~InterfaceTrackingAlg1()
	{}
		
	InterfaceTrackingAlg1::InterfaceTrackingAlg1()
	{
	}
	
	InterfaceTrackingAlg1::InterfaceTrackingAlg1(Data df, DataGen v):
		LevelSetLinear(df),
		velocity(v)
	{		
	}	

		
	void InterfaceTrackingAlg1::initVelocityComputation()
	{
		auto & kk(*kernel);
		unsigned int nv(vectorTemplate->vectors.size());
		unsigned int nd(nD(*vectorTemplate));

		TemplateVE velocityT[nd];
		for(unsigned int i(0); i < nd; ++i)
		{
			velocityT[i].init(*velocity, *vectorTemplate, i);
			kk << velocityT[i].initValues;
		}

		lVelocities.resize(nv);
		for(unsigned int i(0); i < nv; ++i)
		{
			copy(subVE(velocityT[0].values,i,i), lVelocities[i]);
			for(unsigned int j(1); j < nd; ++j)
				copy(cat(lVelocities[i],
				         subVE(velocityT[j].values,i,i)),
				     lVelocities[i]);
		}

	}
		
} // asl

