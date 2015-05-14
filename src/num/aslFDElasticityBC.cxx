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


#include "aslFDElasticityBC.h"
#include "aslFDElasticity.h"
#include <data/aslDataWithGhostNodes.h>
#include "acl/aclMath/aclVectorOfElements.h"
#include "acl/acl.h"
#include "acl/aclHardware.h"
#include "acl/aclGenerators.h"
#include "acl/Kernels/aclExpressionContainer.h"
#include "acl/Kernels/aclKernel.h"
#include "acl/Kernels/aclKernelConfigurationTemplates.h"

#include <math/aslTemplateVE.h>
#include <math/aslTemplateVEExtras.h>
#include <aslGenerators.h>
#include <acl/aclMath/aclMathAlg.h>

#include "aslBasicBC.h"


using acl::generateVEConstant;
using acl::generateVEConstantN;

namespace asl
{
	BCRigidWall::BCRigidWall(SPFDElasticityIncompressibleStatic nm):
		BCond(nm->getDisplacementData()->getBlock(),nm->vectorTemplate),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm)
	{
	}


	void BCRigidWall::init()
	{
		unsigned int nC(num->getDisplacementData()->getEContainer().size());
		loadIndicesToACL();
		loadNeighbourIndicesToACL();
		(*kernel)<<(acl::excerpt(num->getDisplacementData()->getEContainer(),
		                         *indicesACL) =
		            acl::generateVEConstantN(nC,0.));
		(*kernel)<<(acl::excerpt(num->getPressureData()->getEContainer(),
		                         *indicesACL) =
		         acl::excerpt(num->getPressureData()->getEContainer(),
		                      *neighbourIndicesACL));
		kernel->setup();
	}


	void BCRigidWall::execute()
	{
		kernel->compute();
	}

	SPBCond generateBCRigidWall(SPFDElasticityIncompressibleStatic nm, 
                                 const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCRigidWall>(nm));
		addSlices(*bc,sl);	
		return bc;
	}

	BCRigidWallRelaxation::BCRigidWallRelaxation(SPFDElasticityRelaxation nm):
		BCond(nm->getDisplacementData()->getBlock(),nm->vectorTemplate),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		value(generateVEConstantN(nD(*nm),0.))
	{
	}

	BCRigidWallRelaxation::BCRigidWallRelaxation(SPFDElasticityRelaxation nm, acl::VectorOfElements v):
		BCond(nm->getDisplacementData()->getBlock(),nm->vectorTemplate),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		value(v)
	{
	}
		

	void BCRigidWallRelaxation::init()
	{
		loadIndicesToACL();
		loadNeighbourIndicesToACL();
		(*kernel)<<(acl::excerpt(num->getDisplacementData()->getEContainer(),
		                         *indicesACL) = value);
		(*kernel)<<(acl::excerpt(num->getPressureData()->getEContainer(),
		                         *indicesACL) =
		         acl::excerpt(num->getPressureData()->getEContainer(),
		                      *neighbourIndicesACL));
		kernel->setup();
	}


	void BCRigidWallRelaxation::execute()
	{
		kernel->compute();
	}

	SPBCond generateBCRigidWall(SPFDElasticityRelaxation nm, 
                                 const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCRigidWallRelaxation>(nm));
		addSlices(*bc,sl);	
		return bc;
	}

	SPBCond generateBCRigidWall(SPFDElasticityRelaxation nm, 
	                             const AVec<> & u0,
                                 const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCRigidWallRelaxation>(nm,acl::generateVEConstant(u0)));
		addSlices(*bc,sl);	
		return bc;
	}
		
	SPBCond generateBCRigidWall(SPFDElasticity2 nm, 
                                 const std::vector<SlicesNames> & sl)
	{
		auto d(nm->getDisplacementData());
		return generateBCConstantValue(d, AVec<>(nD(d->getBlock()),0), sl);		
	}

	SPBCond generateBCRigidWall(SPFDElasticity2 nm, 
	                             const AVec<> & u0,
                                 const std::vector<SlicesNames> & sl)
	{
		auto d(nm->getDisplacementData());
		return generateBCConstantValue(d, u0, sl);		
	}
		
	BCFreeSurface::BCFreeSurface(FDElasticityIncompressibleStatic* nm):
		BCond(nm->getDisplacementData()->getBlock(),nm->vectorTemplate),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm)
	{
	}


	void BCFreeSurface::init()
	{
		loadIndicesToACL();
		loadNeighbourIndicesToACL();
		(*kernel)<<(acl::excerpt(num->getDisplacementData()->getEContainer(),
		                         *indicesACL) =
		         acl::excerpt(num->getDisplacementData()->getEContainer(),
		                      *neighbourIndicesACL));
		kernel->setup();
	}


	void BCFreeSurface::execute()
	{
		kernel->compute();
	}

	BCFreeSurface2::BCFreeSurface2(SPFDElasticity2 nm):
		BCond(nm->getDisplacementData()->getBlock(),nm->vectorTemplate),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm)
	{
	}

	BCFreeSurface2::~BCFreeSurface2()
	{
	}


	void BCFreeSurface2::init()
	{
		loadIndicesToACL();
		loadNeighbourIndicesToACL();
		(*kernel)<<(acl::excerpt(num->getDisplacementData()->getEContainer(),
		                         *indicesACL) =
		         acl::excerpt(num->getDisplacementData()->getEContainer(),
		                      *neighbourIndicesACL));
		kernel->setup();
	}


	void BCFreeSurface2::execute()
	{
		kernel->compute();
	}

		
	BCImposedDisplacementVelocityValue::
	BCImposedDisplacementVelocityValue(SPFDElasticityIncompressibleStatic nm):
		BCond(nm->getDisplacementData()->getBlock(),nm->vectorTemplate),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		displacement(AVec<double>(nD(*nm),0)),
		bDisplacement(false),	
		velocity(AVec<double>(nD(*nm),0)),
		bVelocity(false),	
		initialized(false)
	{
	}

	void BCImposedDisplacementVelocityValue::setDisplacement(AVec<double> d)
	{
		if (initialized && !bDisplacement)
			errorMessage ("Error (BCImposedDisplacementVelocityValue): An attempt to set displacment value after initialization");
		displacement=d;
		bDisplacement=true;
	}


	void BCImposedDisplacementVelocityValue::setVelocity(AVec<double> v)
	{
		if (initialized && !bVelocity)
			errorMessage ("Error (BCImposedDisplacementVelocityValue): An attempt to set velocity value after initialization");
		velocity=v;
		bVelocity=true;
	}
		
	void BCImposedDisplacementVelocityValue::init()
	{
		loadIndicesToACL();
		if (bDisplacement)
			(*kernel)<<(
			            acl::excerpt(num->getDisplacementData()->getEContainer(),*indicesACL)=
			            acl::generateVEVariableR(displacement)
			            );
		if (bVelocity)
//			(*kernel)<<(
//			            acl::excerpt(num->getVelocityData()->getEContainer(),*indicesACL)=
//			            acl::generateVEVariableR(velocity)
//			            );
		kernel->setup();
		initialized=true;
	}


	void BCImposedDisplacementVelocityValue::execute()
	{
		kernel->compute();
	}


	BCAccelerationSource2::BCAccelerationSource2(FDElasticity2* nm):
		BCond(nm->getDisplacementData()->getBlock(),nm->vectorTemplate),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		acceleration(nD(*nm)),
		initialized(false)
	{
	}


	void BCAccelerationSource2::setAcceleration(AVec<double> a)
	{
		copy(acl::generateVEVariableR(a)*(num->getDeltat()*num->getDeltat()),acceleration);
	}
		
	void BCAccelerationSource2::init()
	{
		loadIndicesToACL();
			(*kernel)<<(
			            acl::excerpt(num->getDisplacementInternalData()->getEContainer(),*indicesACL)-=
			            acceleration);
		kernel->setup();
		initialized=true;
	}


	void BCAccelerationSource2::execute()
	{
		kernel->compute();
	}
		
	BCZeroStressMap::BCZeroStressMap(SPAbstractDataWithGhostNodes d, 
			                         acl::VectorOfElements l,
			                         acl::VectorOfElements m,
			                         SPAbstractDataWithGhostNodes map,
			                         const VectorTemplate *const t):
		BCondWithMap(map,t),
		displacement(d),
		lambda(l),
		mu(m),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)) //< Important _BASIC has better performance
		kernel(new acl::Kernel(acl::KERNEL_BASIC)) //< It is important to fix generateDCFullSafe- SIMDUA issue
	{	
	}

	BCZeroStressMap::~BCZeroStressMap()
	{
	}


	void BCZeroStressMap::init()
	{		
		acl::ExpressionContainer kk;

		auto dispX(generateDCFullSafe(displacement));
		initMapInfrastructure(kk);

		auto & vto(*mapTVE->vto);
		unsigned int nd(nD(*templ));		
        vector<TemplateVE> duTVE(nd);
		for(unsigned int i(0); i < nd; ++i)
		{
			duTVE[i].init(*dispX, *templ, i);
			kk << duTVE[i].initValues;
		}
		
		unsigned int nCells(mapTVE->vto->elementaryCells.size()); 
		acl::TypeID type(getElementType(mapTVE->values));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		auto vNormal(generateVEPrivateVariable(nd,type));
		auto vE(generateVEConstantN(nd,1.));
		auto vB(generateVEPrivateVariable(nd,type));
		auto mA(acl::generateMEPrivateVariable(nd,nd,type));
		auto mDelta(acl::generateMEUnit(nd));
		acl::MatrixOfElements mU(nd,nd);
		acl::VectorOfElements u0(nd);
		for(unsigned int j(0); j < nd; ++j)
			u0[j]=duTVE[j].values[0];
		auto muDlambda(mu/lambda); 
		auto vU0(generateVEPrivateVariable(nd,type));
		auto vU0Res(generateVEPrivateVariable(nd,type));
		auto isMatrixSol(generateVEPrivateVariable(1,typeI));
		auto isBoundary(generateVEPrivateVariable(1,typeI));

		kk << (vNormal = generateVEConstantN(nd, 0.));
		kk << (mA = generateMEDiagonal(generateVEConstantN(nd, 0.)));
		kk << (vB = generateVEConstantN(nd, 0.));
		for(unsigned int i(0); i < nCells; ++i)
		{
			auto cellVal(cellValues(*mapTVE, i));
			kk << (isBoundary = (isGhostNode(0))  && 
			                    (sumOfElements(sign(subVE(cellVal,1,nd))) > 
			                     ((double)nd - .2)));
			kk << (vNormal += select(generateVEConstantN(nd,0.),
			                         gradient(*mapTVE, i),
			                         isBoundary, 
			                         type));
		}
		kk << (gcNormalize(vNormal));
		
		for(unsigned int i(0); i < nCells; ++i)
		{
			auto mT(acl::generateMEConstant(vto.cellMatrices[i]));

			for(unsigned int j(0); j < nd; ++j)
				for(unsigned int k(0); k < nd; ++k)
					mU.setElement(j,k, subVE(duTVE[j].values, 
					                         vto.elementaryCells[i][k+1])[0]);

			auto cellVal(cellValues(*mapTVE, i));
			kk << (isBoundary = (isGhostNode(0))  && 
			                    (sumOfElements(sign(subVE(cellVal,1,nd))) > 
			                     ((double)nd - .2)) ); 

			kk << (mA.getInternalVector() += select((elementProduct(vNormal, vE * mT) + 
			                                         elementProduct(vE * mT,  muDlambda * vNormal) + 
			                                         mDelta * generateME( vE * mT * (muDlambda * vNormal))).getInternalVector(),
			                                        isBoundary,
			                                        type));
			kk << (vB += select(trace(mU, mT) * vNormal +
			                    (muDlambda * vNormal * mU) * mT + 
			                    mU * (mT * vNormal* muDlambda),
			                    isBoundary,
			                    type));
		}
		kk << gcSolveSystem(mA,vB,vU0);
		kk << (isMatrixSol = (acl::fabs(det(mA)) > .05));
		
		auto counterD(generateVEPrivateVariable(1,type));
		kk << (vU0 = generateVEConstantN(nd, 0.));
		kk << (counterD = generateVEConstant(0));
		unsigned int nDir(mapTVE->vectorTemplate->vectors.size());
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = (isComputationNode(i)));
			kk << (counterD += select(generateVEConstant(1.), isBoundary, type));
			acl::VectorOfElements vv(nd); 
			for(unsigned int k(0); k < nd; ++k)
				vv[k] = duTVE[k].values[i];
			kk << ( vU0 += select(vv, isBoundary, type));			
		}
		kk << (vU0 /= max(counterD, generateVEConstant(.1), type));
		kk << (vU0Res = select(vU0, vU0Res, isMatrixSol, type));
		kk << (vU0Res = select(u0, vU0Res, isGhostNode(0), type));
		
		kk << (acl::assignmentSafe(displacement->getEContainer(), vU0Res));

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((isGhostNode() && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));
		kernel->setup();
//		cout<<kernel->getKernelSource()<<endl;

	}
/*
	void BCZeroStressMap::init()
	{		
		acl::ExpressionContainer kk;

		auto dispX(generateDCFullSafe(displacement));
		initMapInfrastructure(kk);

		auto & vto(*mapTVE->vto);
		unsigned int nd(nD(*templ));		
        vector<TemplateVE> duTVE(nd);
		for(unsigned int i(0); i < nd; ++i)
		{
			duTVE[i].init(*dispX, *templ, i);
			kk << duTVE[i].initValues;
		}
		
		unsigned int nCells(mapTVE->vto->elementaryCells.size()); 
		acl::TypeID type(getElementType(mapTVE->values));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		auto vNormal(generateVEPrivateVariable(nd,type));
		auto vE(generateVEConstantN(nd,1.));
		auto vB(generateVEPrivateVariable(nd,type));
		auto mA(acl::generateMEPrivateVariable(nd,nd,type));
		auto mDelta(acl::generateMEUnit(nd));
		acl::MatrixOfElements mU(nd,nd);
		acl::VectorOfElements u0(nd);
		for(unsigned int j(0); j < nd; ++j)
			u0[j]=duTVE[j].values[0];
		auto muDlambda(mu/lambda); 
		auto vU0(generateVEPrivateVariable(nd,type));
		auto vU0Res(generateVEPrivateVariable(nd,type));
		auto counter(generateVEPrivateVariable(1,type));
		auto isBoundary(generateVEPrivateVariable(1,typeI));

		kk << (counter = generateVEConstant(0));
		kk << (vU0Res = generateVEConstantN(nd, 0.));
		kk << (vNormal = generateVEConstantN(nd, 0.));
		for(unsigned int i(0); i < nCells; ++i)
		{
			auto cellVal(cellValues(*mapTVE, i));
			kk << (isBoundary = (isGhostNode(0))  && 
			                    (sumOfElements(sign(subVE(cellVal,1,nd))) > 
			                     ((double)nd - .2)));
			kk << (vNormal += select(generateVEConstantN(nd,0.),
			                         gradient(*mapTVE, i),
			                         isBoundary, 
			                         type));
		}
		kk << (gcNormalize(vNormal));
		
		for(unsigned int i(0); i < nCells; ++i)
		{
			auto mT(acl::generateMEConstant(vto.cellMatrices[i]));

			for(unsigned int j(0); j < nd; ++j)
				for(unsigned int k(0); k < nd; ++k)
					mU.setElement(j,k, subVE(duTVE[j].values, 
					                         vto.elementaryCells[i][k+1])[0]);
//			kk << (vNormal = gradient(*mapTVE, i));
//			kk << (gcNormalize(vNormal));
			kk << (mA = elementProduct(vNormal, vE * mT) + 
			            elementProduct(vE * mT,  muDlambda * vNormal) + 
			            mDelta * generateME( vE * mT * (muDlambda * vNormal)));
			kk << (vB =  trace(mU, mT) * vNormal +
			             (muDlambda * vNormal * mU) * mT + 
			            mU * (mT * vNormal* muDlambda));
			
			kk << gcSolveSystem(mA,vB,vU0);
			auto cellVal(cellValues(*mapTVE, i));
			kk << (isBoundary = (isGhostNode(0))  && 
			                    (sumOfElements(sign(subVE(cellVal,1,nd))) > 
			                     ((double)nd - .2)) && (acl::fabs(det(mA)) > .1));
			kk << (counter += select(generateVEConstant(0.),
			                        generateVEConstant(1.),
			                        isBoundary, 
			                        type));
			kk << ( vU0Res += select(generateVEConstantN(nd,0.),
			                         vU0,
			                         isBoundary, 
			                         type));
		}
		kk << (vU0Res /= max(counter, generateVEConstant(.1), type));

		auto counterD(generateVEPrivateVariable(1,type));
		kk << (vU0 = generateVEConstantN(nd, 0.));
		kk << (counterD = generateVEConstant(0));
		unsigned int nDir(mapTVE->vectorTemplate->vectors.size());
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = (isComputationNode(i)));
			kk << (counterD += select(generateVEConstant(0.),
			                        generateVEConstant(1.),
			                        isBoundary, 
			                        type));
			acl::VectorOfElements vv(nd); 
			for(unsigned int k(0); k < nd; ++k)
				vv[k] = duTVE[k].values[i];
			kk << ( vU0 += select(generateVEConstantN(nd,0.),
			                      vv,
			                      isBoundary, 
			                      type));			
		}
		kk << (vU0 /= max(counterD, generateVEConstant(.1), type));
		kk << (vU0Res = select(vU0, vU0Res, counter > .1, type));
		kk << (vU0Res = select(u0,
		                       vU0Res,
		                       isGhostNode(0),
		                       type));
		
		kk << (acl::assignmentSafe(displacement->getEContainer(), vU0Res));

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((isGhostNode() && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));
		kernel->setup();
//		cout<<kernel->getKernelSource()<<endl;

	}
*/

	void BCZeroStressMap::execute()
	{
		kernel->compute();
	}		
		
	SPNumMethod generateBCZeroStress(SPElasticityCommonA nm, 
                                     SPAbstractDataWithGhostNodes map)
	{
		return std::make_shared<BCZeroStressMap>(
		                       nm->getDisplacementData(),
		                       nm->getBulkModulus() - nm->getShearModulus()/1.5,
		                       nm->getShearModulus(),
		                       map,
		                       nm->vectorTemplate);
	}

	SPNumMethod generateBCZeroStressP(SPFDElasticityIncompressibleStatic nm, 
                                      SPAbstractDataWithGhostNodes map)
	{
		return generateBCConstantValueMiddlePoint(nm->getPressureData(), 
	                                 			  0, 
	                                  			  map,
	                                  			  nm->vectorTemplate);
	}

	SPNumMethod generateBCZeroStressP(SPFDElasticityRelaxation nm, 
                                      SPAbstractDataWithGhostNodes map)
	{
		return generateBCConstantValueMiddlePoint(nm->getPressureData(), 
	                                 			  0, 
	                                  			  map,
	                                  			  nm->vectorTemplate);
	}

		
} // asl

