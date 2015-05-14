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


#include "aslFDPoroElasticityBC.h"
#include "aslFDPoroElasticity.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include "aslFDElasticityBC.h"
#include "aslBasicBC.h"
#include "acl/Kernels/aclKernelConfigurationTemplates.h"
#include "acl/Kernels/aclKernel.h"
#include <acl/acl.h>
#include <acl/aclGenerators.h>
#include <acl/Kernels/aclExpressionContainer.h>
#include <math/aslTemplateVE.h>
#include <aslGenerators.h>
#include <math/aslIndex2Position.h>

namespace asl
{

	BCRigidWallPoroElasticity::BCRigidWallPoroElasticity(SPFDPoroElasticity nm):
		BCond(nm->getDisplacementData()->getBlock(),nm->vectorTemplate),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		value(acl::generateVEConstantN(nD(*nm),0.))
	{
	}

	BCRigidWallPoroElasticity::~BCRigidWallPoroElasticity(){}

	BCRigidWallPoroElasticity::BCRigidWallPoroElasticity(SPFDPoroElasticity nm, acl::VectorOfElements v):
		BCond(nm->getDisplacementData()->getBlock(),nm->vectorTemplate),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		value(v)
	{
	}
		

	void BCRigidWallPoroElasticity::init()
	{
		loadIndicesToACL();
		loadNeighbourIndicesToACL();
		(*kernel)<<(acl::excerpt(num->getDisplacementData()->getEContainer(),
		                         *indicesACL) = value);
		(*kernel)<<(acl::excerpt(num->getPressureData()->getEContainer(),
		                         *indicesACL) =
		         acl::excerpt(num->getPressureData()->getEContainer(),
		                      *neighbourIndicesACL));
		(*kernel)<<(acl::excerpt(num->getLiquidPressureData()->getEContainer(),
		                         *indicesACL) =
		            acl::excerpt(num->getLiquidPressureData()->getEContainer(),
		                         *neighbourIndicesACL));
		kernel->setup();
	}


	void BCRigidWallPoroElasticity::execute()
	{
		kernel->compute();
	}

	BCRigidWallDF::BCRigidWallDF(SPFDPoroElasticity nm, 
	                             SPDistanceFunction rw, 
	                             SPAbstractDataWithGhostNodes map):
		BCondWithMap(map,nm->vectorTemplate),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)), //< It is important to fix generateDCFullSafe- SIMDUA issue
		num(nm),
		rWall(rw)
	{	
	}

	BCRigidWallDF::~BCRigidWallDF()
	{
	}


	void BCRigidWallDF::init()
	{		
		double wallVal(.0);
		acl::ExpressionContainer kk;
		initMapInfrastructure(kk);

		auto dispX(generateDCFullSafe(num->getDisplacementData()));
		auto pressureX(generateDCFullSafe(num->getPressureData()));
		auto pressureLX(generateDCFullSafe(num->getLiquidPressureData()));

		acl::TypeID type(getElementType(mapTVE->values));
		
//		unsigned int nd(nD(*templ));
		unsigned int nDir(mapTVE->values.size());
//		double dx(map->getBlock().dx);
		auto counter(generateVEPrivateVariable(1,type));
		auto val(generateVEPrivateVariable(1,type));
		auto inWall(generateVEPrivateVariable(1,type));
		Index2PositionACL i2p(map->getBlock(), type);
			
		kk <<(inWall=rWall->getDistance(i2p.positionWithInit +
		                                num->getDisplacementData()->getEContainer()));
		kk << (counter = acl::generateVEConstant(0));

		TemplateVE duTVE;
		duTVE.init(*pressureX, *templ, 0, false);
		kk << duTVE.initValues;
		kk << (val = acl::generateVEConstant(0));
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (counter += select(acl::generateVEConstant(1.), isComputationNode(i), type));
			kk << (val += select(duTVE.getValue(i), isComputationNode(i), type));			
		}
		kk << (counter = max(counter, acl::generateVEConstant(.1), type));
		kk << (acl::assignmentSafe(num->getPressureData()->getEContainer(), 
		                           select(num->getPressureData()->getEContainer(),
		                                  val/counter,
		                                  inWall > wallVal,
		                                  type)));

		duTVE.init(*pressureLX, *templ, 0, false);
		kk << duTVE.initValues;
		kk << (val = acl::generateVEConstant(0));
		for(unsigned int i(1); i < nDir; ++i)
			kk << (val += select(duTVE.getValue(i), isComputationNode(i), type));			
		kk << (acl::assignmentSafe(num->getLiquidPressureData()->getEContainer(), 
		                           select(num->getLiquidPressureData()->getEContainer(),
		                                  val/counter,
		                                  inWall > wallVal,
		                                  type)));
/*
		for(unsigned int i(0); i < nd; ++i)
		{
			duTVE.init(*dispX, *templ, i, false);
			kk << (val = acl::generateVEConstant(0));		
			kk << duTVE.initValues;
			for(unsigned int j(1); j < nDir; ++j)
				kk << (val += select(duTVE.getValue(j)+(.5*dx*templ->vectors[j][i]*(inWall-wallVal)), 
				                     isComputationNode(j), type));			
			kk << (acl::assignmentSafe(subVE(num->getDisplacementData()->getEContainer(),i), 
			                           select(subVE(num->getDisplacementData()->getEContainer(),i),
			                                  val/counter,
			                                  inWall > wallVal,
			                                  type)));
		}
*/
		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((isGhostNode() && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));
		kernel->setup();
	}

	void BCRigidWallDF::execute()
	{
		kernel->compute();
	}		
		
	void addBCRigidWall(std::vector<SPNumMethod> & bcList,
	                     SPFDPoroElasticity nm, 
	                     const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCRigidWallPoroElasticity>(nm));
		addSlices(*bc,sl);	
		bcList.push_back(bc);		
	}

	void addBCZeroStress(std::vector<SPNumMethod> & bcList,
	                     SPFDPoroElasticity nm, 
	                     SPAbstractDataWithGhostNodes map)
	{
		bcList.push_back(asl::generateBCZeroStress(nm, map));
		bcList.push_back(generateBCConstantValueMiddlePoint(nm->getPressureData(), 
		                                                    0, 
		                                                    map,
		                                                    nm->vectorTemplate));
		bcList.push_back(generateBCConstantValueMiddlePoint(nm->getLiquidPressureData(), 
		                                                    0, 
		                                                    map,
		                                                    nm->vectorTemplate));
	}

	void addBCZeroStress(std::vector<SPNumMethod> & bcList,
	                     SPFDPoroElasticity nm,
	                     SPPositionFunction p,
	                     SPAbstractDataWithGhostNodes map)
	{
		bcList.push_back(asl::generateBCZeroStress(nm, map));
		bcList.push_back(generateBCConstantValue(nm->getPressureData(), 0, map));
		bcList.push_back(generateBCConstantValue(nm->getLiquidPressureData(), p, map));
	}
		
	void addBCRigidWallDF(std::vector<SPNumMethod> & bcList,
	                      SPFDPoroElasticity nm,
	                      SPDistanceFunction rw, 
	                      SPAbstractDataWithGhostNodes map)
	{
		auto bc(std::make_shared<BCRigidWallDF>(nm, rw, map)); 
		bcList.push_back(bc);
	}


	void addBCRigidWallDF(std::vector<SPNumMethod> & bcList,
	                      SPFDPoroElasticity nm,
	                      SPAbstractDataWithGhostNodes rw, 
			              SPAbstractDataWithGhostNodes map)
	{
		auto rwX(generateDCFullSafe(rw,-1.));
		auto rwXX(std::make_shared<DataInterpolation>(rwX));
		addBCRigidWallDF(bcList, nm, rwXX, map);		
	}
} // asl

