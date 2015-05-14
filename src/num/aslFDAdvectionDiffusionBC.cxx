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


#include "aslFDAdvectionDiffusionBC.h"
#include "aslFDAdvectionDiffusion.h"
#include <aslGenerators.h>
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include <acl/aclHardware.h>
#include <math/aslTemplateVE.h>

#include <data/aslDataWithGhostNodes.h>
#include <acl/Kernels/aclKernel.h>
#include <math/aslTemplates.h>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>
#include <math/aslDistanceFunctionAlg.h>

namespace asl
{

	BCConstantFluxMap::BCConstantFluxMap(Data d,
	                                     const acl::VectorOfElements & val, 
	                                     Data map,
	                                     const VectorTemplate *const t):
		BCondWithMap(map, t),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		data(d),
		value(val)
	{
	}

	BCConstantFluxMap::~BCConstantFluxMap()
	{
	}


	void BCConstantFluxMap::init()
	{		
		acl::ExpressionContainer kk;

		unsigned int nd(nD(data->getBlock())); //< number of dimmentions
		unsigned int nc(data->getEContainer().size()); // number of components 
		
		auto dataX(generateDCFullSafe(data));
		
		initMapInfrastructure(kk);
					
		acl::TypeID type(getElementType(data->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		TemplateVE dataTVE(*dataX, *templ);
		kk << dataTVE.initValues;
		
		auto normal(generateVEPrivateVariable(nd,type));		
		auto area(generateVEPrivateVariable(1,type));		
		auto sw(generateVEPrivateVariable(1,type));
		auto swc(generateVEPrivateVariable(1,type));
		auto sna(generateVEPrivateVariable(1,type));		
		auto sna1mx(generateVEPrivateVariable(1,type));		
		auto snaxc(generateVEPrivateVariable(nc,type));
		auto isBoundary(generateVEPrivateVariable(1,typeI));
		
		unsigned int nDir(templ->vectors.size());
		
		kk << (sw = acl::generateVEConstant(0.));
		kk << (swc = acl::generateVEConstant(0.));
		kk << (gcNormalize(normal));
//		kk << gcBoundaryAreaPerGhostPoint(*mapTVE, area);

//		TemplateVE mSTVE(mapTVE->values*.5+catN(mapTVE->getValue(0),nDir)*.5 ,
//		                 *mapTVE->vectorTemplate);
//		kk << gcBoundaryArea(mSTVE, area); 
		kk << gcBoundaryArea(*mapTVE, area); 

		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isComputationNode(i));
			kk << (sw += select(acl::generateVEConstant(templ->laplasCoefs[i]), 
			                    isBoundary,type));
			kk << (swc += select(dataTVE.getValue(i) * templ->laplasCoefs[i], 
			                    isBoundary,type));
		}
		auto vLocal((swc+sw*value)/sw);
//		auto vLocal((swc+area*value)/sw);
//		auto vLocal(area*.25);
//		kk << (acl::assignmentSafe(data->getEContainer(), vLocal));
		kk << (acl::assignmentSafe(data->getEContainer(), 
		                           select(dataTVE.getValue(0), vLocal, isGhostNode(0), type)));
//		
/*		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((isGhostNode() && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));
*/  
		(*kernel)<<kk;
		kernel->setup();
	}

	void BCConstantFluxMap::execute()
	{
		kernel->compute();
	}


		
		
	SPNumMethod generateBCConstantFlux(SPFDAdvectionDiffusion nm, 
	                                   double flux,
	                                   SPAbstractDataWithGhostNodes map)
	{
		auto a(make_shared<BCConstantFluxMap>(nm->getData()[0], 
		                                      flux/nm->getDiffusionCoefficient(),
		                                      map,
		                                      nm->getVectorTemplate()));
		return a;
		
	}
		
} //asl
