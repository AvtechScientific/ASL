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


#include "aslBasicBC2.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include "acl/acl.h"
#include "acl/aclGenerators.h"
#include <acl/Kernels/aclKernel.h>
#include "acl/Kernels/aclKernelConfigurationTemplates.h"
#include <aslGenerators.h>
#include <math/aslTemplateVE.h>
#include "acl/aclMath/aclMathAlg.h"
#include <math/aslIndex2Position.h>
#include <math/aslPositionFunction.h>
#include <utilities/aslUValue.h>
#include <math/aslDistanceFunctionAlg.h>


namespace asl
{

/*	BCConstantValueMap::BCConstantValueMap(Data d, 
	                                       const acl::VectorOfElements & v, 
	                                       Data map):
		BCondWithMap(map, nearestNeigboursVT(nD(d->getBlock()))),
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		data(d),
		value(v)
	{
	}

	BCConstantValueMap::~BCConstantValueMap()
	{
	}


	void BCConstantValueMap::init()
	{		
		acl::ExpressionContainer kk;
		
		acl::TypeID type(getElementType(data->getEContainer()));

		auto vGhost(data->getEContainer());
		kk << (acl::assignmentSafe(vGhost, 
		                           select(vGhost, value, isGhostNode(), type)));

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((map->getEContainer() <= 0 && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}

	void BCConstantValueMap::execute()
	{
		kernel->compute();
	}
*/
/*
	BCValuePFMap::BCValuePFMap(Data d, 
	                           SPPositionFunction v, 
	                           Data map):
		BCondWithMap(map, nearestNeigboursVT(nD(d->getBlock()))),
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		data(d),
		value(v)
	{
	}

	BCValuePFMap::~BCValuePFMap()
	{
	}

	void BCValuePFMap::init()
	{		
		acl::ExpressionContainer kk;
		
		acl::TypeID type(getElementType(data->getEContainer()));

		auto vGhost(data->getEContainer());
		asl::Index2PositionACL i2p(data->getBlock(),type);

		kk << i2p.initPosition;		
		auto val(value->value(i2p.position));
		kk << (acl::assignmentSafe(vGhost, 
		                           select(vGhost, val, isGhostNode(), type)));

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((map->getEContainer() <= 0 && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}

	void BCValuePFMap::execute()
	{
		kernel->compute();
	}
*/
		
	BCConstantGradientMap2::BCConstantGradientMap2(Data d, 
	                                               const acl::VectorOfElements & v, 
	                                               Data map,
	                                               const VectorTemplate *const t):
		BCondWithMap(map, t),
//		kernelCN(new acl::Kernel(acl::KERNEL_SIMDUA)),
//		kernelGN(new acl::Kernel(acl::KERNEL_SIMDUA)),
		kernelCN(new acl::Kernel(acl::KERNEL_BASIC)),
		kernelGN(new acl::Kernel(acl::KERNEL_BASIC)),
		data(d),
		value(v)
	{
	}

	BCConstantGradientMap2::BCConstantGradientMap2(Data d, 
	                                     const acl::VectorOfElements & v, 
	                                     Data map,
	                                     Data computationalDomain,
	                                     const VectorTemplate *const t):
		BCondWithMap(map, computationalDomain, t),
//		kernelCN(new acl::Kernel(acl::KERNEL_SIMDUA)),
//		kernelGN(new acl::Kernel(acl::KERNEL_SIMDUA)),
		kernelCN(new acl::Kernel(acl::KERNEL_BASIC)),
		kernelGN(new acl::Kernel(acl::KERNEL_BASIC)),
		data(d),
		value(v)
	{
	}
		
	BCConstantGradientMap2::~BCConstantGradientMap2()
	{
	}


	void BCConstantGradientMap2::init()
	{		
		unsigned int nd(nD(data->getBlock()));
		unsigned int nc(data->getEContainer().size()); // number of components 
		unsigned int nDir(templ->vectors.size());
		acl::TypeID type(getElementType(data->getEContainer()));
		
		auto counter(generateVEPrivateVariable(1,type));
		auto vLocal(generateVEPrivateVariable(1,type));
		auto normal(generateVEPrivateVariable(nd,type));		
		auto dirCoef(generateVEPrivateVariable(1,type));
		auto xx(generateVEPrivateVariable(1,type));

		auto dataX(generateDCFullSafe(data));
		TemplateVE dataTVE;	

		acl::ExpressionContainer kkCN;
		initMapInfrastructure(kkCN);
		kkCN << (normal = gradient(*mapTVE));
		kkCN << (gcNormalize(normal));
		kkCN << (counter = acl::generateVEConstant(0.));
		kkCN << (xx = acl::generateVEConstant(0.));

		for(unsigned int i(1); i < nDir; ++i)
		{
			kkCN << (dirCoef = (normal * normalize(templ->vectors[i])) + 1.);
			kkCN << (dirCoef *= dirCoef );			
			kkCN << (counter+= select(dirCoef, isGhostNode(templ->invertVectors[i]), type));
			kkCN << (xx+= select(dirCoef * exBoundaryX(*mapTVE,templ->invertVectors[i]), 
			                     isGhostNode(templ->invertVectors[i]), type));
		}
		kkCN << (xx=select(acl::generateVEConstant(1.),xx/counter,counter>0.001,type));
		kkCN << (counter = max(counter,acl::generateVEConstant(.001), type));
		
		
		for(unsigned int ii(0); ii < nc; ++ii)
		{
			dataTVE.init(*dataX, *templ, ii,false);
			kkCN << (vLocal = acl::generateVEConstant(0.));
			for(unsigned int i(1); i < nDir; ++i)
			{
				kkCN << (dirCoef = (normal * normalize(templ->vectors[i])) + 1.);
				kkCN << (dirCoef *= dirCoef );							
				kkCN << (vLocal += select((dataTVE.getValue(i) - 
				                           (normal*templ->vectors[i])*subVE(value, ii)) *
				                          dirCoef, 
				                          isGhostNode(templ->invertVectors[i]),
				                          type));
			}
			kkCN << (acl::assignmentSafe(subVE(data->getEContainer(),ii), 
			                             select(dataTVE.getValue(0),
			                             xx*dataTVE.getValue(0) + (1.-xx)*vLocal/counter, 
			                                    isComputationNode(0), 
			                                    type)));
		}
		kernelCN->addExpression(acl::elementOperators::
			                        ifElse(acl::elementOperators::any((isComputationNode() && 
			                                                           map->getEContainer() < .9999)[0]),
			                               kkCN.expression, 
			                               {} ));
		kernelCN->setup();

		
		acl::ExpressionContainer kkGN;
		initMapInfrastructure(kkGN);
		kkGN << (normal = gradient(*mapTVE));
		kkGN << (gcNormalize(normal));
		kkGN << (counter = acl::generateVEConstant(0.));
					
		
		for(unsigned int i(1); i < nDir; ++i)
		{
			kkGN << (dirCoef = (normal * normalize(templ->vectors[i])) + 1.);
			kkGN << (dirCoef *= dirCoef );			
			kkGN << (counter+= select(dirCoef, isComputationNode(i), type));
		}
		kkGN << (counter = max(counter,acl::generateVEConstant(.001), type));

		for(unsigned int ii(0); ii < nc; ++ii)
		{
			dataTVE.init(*dataX, *templ, ii,false);
			kkGN << (vLocal = acl::generateVEConstant(0.));
			for(unsigned int i(1); i < nDir; ++i)
			{
				kkGN << (dirCoef = (normal * normalize(templ->vectors[i])) + 1.);
				kkGN << (dirCoef *= dirCoef );							
				kkGN << (vLocal += select((dataTVE.getValue(i) - 
				                           (normal*templ->vectors[i])*subVE(value, ii)) *
				                          dirCoef, 
				                          isComputationNode(i),
				                          type));
			}
			kkGN << (acl::assignmentSafe(subVE(data->getEContainer(),ii), 
			                             select(subVE(dataTVE.values,0), 
			                                    vLocal/counter, 
			                                    isGhostNode(0), 
			                                    type)));
		}
		kernelGN->addExpression(acl::elementOperators::
			                        ifElse(acl::elementOperators::any((isGhostNode() && 
			                                                           map->getEContainer() > -.9999)[0]),
			                               kkGN.expression, 
			                               {} ));
		kernelGN->setup();
	}

	void BCConstantGradientMap2::execute()
	{
		kernelCN->compute();
		kernelGN->compute();
	}

		
} // asl

