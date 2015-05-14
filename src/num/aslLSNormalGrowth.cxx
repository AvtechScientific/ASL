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


#include "aslLSNormalGrowth.h"
#include <math/aslTemplates.h>
#include <data/aslDataWithGhostNodes.h>
#include "acl/acl.h"
#include "acl/Kernels/aclKernel.h"
#include "acl/Kernels/aclKernelConfigurationTemplates.h"
#include "acl/Kernels/aclKernelMerger.h"
#include "acl/aclHardware.h"
#include "acl/aclGenerators.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include "acl/aclMath/aclBarycentric.h"
#include <algorithm>
#include <math/aslTemplates.h>
#include <math/aslTemplateVE.h>
#include <math/aslTemplatesExtras.h>
#include <math/aslTemplateVEExtras.h>
#include <math/aslDistanceFunction.h>


namespace asl
{
	using acl::generateVEConstant;
	using acl::generateVEConstantN;

	LSNormalGrowth::~LSNormalGrowth()
	{}
		
	LSNormalGrowth::LSNormalGrowth()
	{
	}
	
	LSNormalGrowth::LSNormalGrowth(Data df, DataGen c):
		LevelSetLinear(df),
		superSaturation(c)
	{		
	}	

	inline bool isIn(int i, AVec<int> a)
	{
		bool b(false);
		for(unsigned int j(0); j < a.getSize(); ++j )
			b|=(a[j]==i);
		return b;
	}
		
	void LSNormalGrowth::initVelocityComputation()
	{
		auto & kk(*kernel);
		acl::TypeID type(getElementType(distanceField->getEContainer()));
		unsigned int nd(nD(*vectorTemplate));		
		unsigned int nv(vectorTemplate->vectors.size());
		auto vto(distanceTVE->vto);
		unsigned int nCells(vto->elementaryCells.size());


		TemplateVE sSatT(*superSaturation, *vectorTemplate);
		kk << sSatT.initValues;

		lVelocities.resize(nv);

		auto isBoundary(generateVEPrivateVariable(1,acl::TYPE_SELECT[type]));
		auto normal(generateVEPrivateVariable(nd,type));
//		auto invlen(generateVEPrivateVariable(1,type));
		auto counter(generateVEPrivateVariable(1, type));

		kk << (counter = generateVEConstant(0.));
		
		for(unsigned int i(0); i < nv; ++i)
		{
			copy(generateVEPrivateVariable(nd,type), lVelocities[i]);
			kk << (lVelocities[i] = generateVEConstantN(nd,0.));
		}

/*		vector<acl::VectorOfElements> normalC;
		kk << gcGradientAllCells(*distanceTVE,  normalC);
		for(unsigned int j(0); j < nCells; ++j)
			kk << gcNormalize(normalC[j]);
		
		for(unsigned int i(1); i < nv; ++i)
		{
			kk << (normal = generateVEConstantN(nd,0.));
			for(unsigned int j(0); j < nCells; ++j)
				if(isIn(i, vto->elementaryCells[j]))
					kk << (normal += normalC[j]);
			kk << gcNormalize(normal);
			
			kk << (isBoundary = isBoundaryDir(i));
			
			kk << (lVelocities[i] = select(generateVEConstantN(nd,0.), 
			                               normal * getValueOnBoundary(sSatT.values, i),
			                               isBoundary,
			                               type));
			kk << (counter += select(generateVEConstant(0.), 
			                         generateVEConstant(1.),
			                         isBoundary,
			                         type)); 
			kk << (lVelocities[0] += lVelocities[i]);
		}
		kk << (counter = max(counter, generateVEConstant(1.), type));
		kk << (lVelocities[0] = lVelocities[0] / counter);*/
	}
		
} // asl

