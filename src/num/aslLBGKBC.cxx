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


#include "aslLBGKBC.h"

#include "aslLBGK.h"
#include "aslBasicBC.h"
#include <data/aslDataWithGhostNodes.h>
#include "acl/acl.h"
#include "acl/Kernels/aclKernel.h"
#include "acl/Kernels/aclKernelConfigurationTemplates.h"
#include "acl/Kernels/aclKernelMerger.h"
#include "acl/aclGenerators.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include <algorithm>
#include <math/aslTemplates.h>
#include <aslGenerators.h>
#include <math/aslTemplateVE.h>
#include "acl/aclMath/aclMathAlg.h"
#include <math/aslIndex2Position.h>
#include <math/aslPositionFunction.h>

namespace asl
{
	using acl::generateVEConstant;
	using acl::generateVEConstantN;
	using acl::generateVEPrivateVariable;
	

	BCLBGKCommon::BCLBGKCommon(SPLBGK nm):
		BCond(nm->getF()->getBlock(),nm->getVectorTemplate()),
		num(nm),
		kernels(nm->getVectorTemplate()->vectors.size()),
		km(new acl::KernelMerger()),
		directionGroupsShifts(nm->getVectorTemplate()->vectors.size()),
		directionGroupsSizes(nm->getVectorTemplate()->vectors.size())
	{
		for(unsigned int i(0); i<kernels.size(); ++i)
			kernels[i].reset(new acl::Kernel(acl::KERNEL_BASIC));
	}


	void BCLBGKCommon::sortDirections()
	{
		sortTwoVectors(directions, indices);

		unsigned int n(directionGroupsShifts.getSize());
		for(unsigned int i(0); i<n; ++i){
			directionGroupsShifts[i]=std::find(directions.begin(), 
			                                   directions.end(),
			                                   (int)i)-directions.begin();
			directionGroupsSizes[i]=directions.rend()-
				std::find(directions.rbegin(), directions.rend(), (int)i)
				-directionGroupsShifts[i];
			directionGroupsSizes[i]=directionGroupsSizes[i]<0 ? 0 : directionGroupsSizes[i];
		}
	}


	void BCLBGKCommon::execute()
	{
		km->compute();
	}

		
	BCNoSlip::BCNoSlip(SPLBGK nm):
		BCLBGKCommon(nm)
	{}


	void BCNoSlip::init()
	{

		unsigned int nC(directionGroupsShifts.getSize());		
		sortDirections();
		loadIndicesToACL();
		
		for(unsigned int i(1); i<nC; ++i ){
			if(directionGroupsSizes[i]>0)
			{
				auto vGhost(subVE(num->getF()->getEContainer(), i));
				auto vBulk(subVE(num->getF()->getEContainer(), templ->invertVectors[i]));
				int directionShift(block.c2i(templ->vectors[i]));
				acl::VectorOfElements internalInd(acl::generateVEIndex(directionGroupsSizes[i])+
				                                  directionGroupsShifts[i]);
				(*(kernels[i]))<<(acl::excerpt(acl::excerpt(vGhost,*indicesACL)=
					                          acl::excerpt(vBulk,*indicesACL+directionShift),
				                              internalInd));
				km->addKernel(kernels[i]);                
			}
		}
		km->setup();
	}


	BCConstantPressure::BCConstantPressure(SPLBGK nm, const acl::VectorOfElements & p):
		BCLBGKCommon(nm),
		pressure(p)
	{}
		
	void BCConstantPressure::init()
	{

		unsigned int nC(directionGroupsShifts.getSize());		
		sortDirections();
		loadIndicesToACL();
		
		for(unsigned int i(1); i<nC; ++i ){
			if(directionGroupsSizes[i]>0)
			{
				acl::VectorOfElements vGhost(1);
				vGhost[0]=num->getF()->getEContainer()[i];
				auto internalInd(acl::generateVEIndex(directionGroupsSizes[i])+
				                                  directionGroupsShifts[i]);
				(*(kernels[i]))<<(acl::excerpt(acl::excerpt(vGhost,*indicesACL)=
					                          pressure,
				                              internalInd));
				km->addKernel(kernels[i]);                
			}
		}
		km->setup();
	}

	BCConstantVelocity::BCConstantVelocity(SPLBGK nm, const acl::VectorOfElements & v):
		BCLBGKCommon(nm),
		velocity(v)
	{}
		
	void BCConstantVelocity::init()
	{

		unsigned int nC(directionGroupsShifts.getSize());		
		sortDirections();
		loadIndicesToACL();
		
		for(unsigned int i(1); i<nC; ++i ){
			if(directionGroupsSizes[i]>0)
			{
				auto vGhost(subVE(num->getF()->getEContainer(), i));
				auto vBulk(subVE(num->getF()->getEContainer(), templ->invertVectors[i]));
				int directionShift(block.c2i(templ->vectors[i]));
				acl::VectorOfElements internalInd(acl::generateVEIndex(directionGroupsSizes[i])+
				                                  directionGroupsShifts[i]);
				(*(kernels[i]))<<(acl::excerpt(acl::excerpt(vGhost,*indicesACL) =
					                           acl::excerpt(vBulk,*indicesACL+directionShift) +
				                               6. * (velocity * templ->vectors[i]),
				                               internalInd));
				km->addKernel(kernels[i]);                
			}
		}
		km->setup();
	}

	BCConstantPressureVelocity::BCConstantPressureVelocity(SPLBGK nm,
	                                                       const acl::VectorOfElements & p,
	                                                       const acl::VectorOfElements & v):
		BCLBGKCommon(nm),
		pressure(p),
		velocity(v)
	{}
		
	void BCConstantPressureVelocity::init()
	{

		unsigned int nC(directionGroupsShifts.getSize());		
		sortDirections();
		loadIndicesToACL();
		
		for(unsigned int i(1); i<nC; ++i ){
			if(directionGroupsSizes[i]>0)
			{
				auto vGhost(subVE(num->getF()->getEContainer(), i));
				acl::VectorOfElements internalInd(acl::generateVEIndex(directionGroupsSizes[i])+
				                                  directionGroupsShifts[i]);
				(*(kernels[i]))<<(acl::excerpt(acl::excerpt(vGhost,*indicesACL) =
					                           pressure + 3. * (velocity * templ->vectors[i]),
				                               internalInd));
				km->addKernel(kernels[i]);                
			}
		}
		km->setup();
	}
		
		
	BCNoSlipMap::BCNoSlipMap(SPLBGK nm, 
	                         SPAbstractDataWithGhostNodes map):
		BCondWithMap(map,nm->vectorTemplate),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),//< Important _BASIC has better performance
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm)
	{	
	}

	BCNoSlipMap::~BCNoSlipMap()
	{
	}


	void BCNoSlipMap::init()
	{		
		acl::ExpressionContainer kk;

		auto fX(generateDCFullSafe(num->getF()));

		initMapInfrastructure(kk);

		acl::TypeID type(getElementType(num->getF()->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		auto isBoundary(generateVEPrivateVariable(1,typeI));

		unsigned int nDir(mapTVE->vectorTemplate->vectors.size());
		auto & block(fX->getBlock());

		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isGhostNode(0) && isComputationNode(i));

			auto vGhost(subVE(fX->getEContainer(), i));
			int directionShift(block.c2i(templ->vectors[i]));
			auto vBulk(subVE(fX->getEContainer(), 
			                                             templ->invertVectors[i]));			                                       
			vBulk[0] = acl::generateShiftedElement(vBulk[0],directionShift);
			kk << (acl::assignmentSafe(vGhost, select(vGhost, vBulk, isBoundary, type)));
		}		

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((map->getEContainer() <= 0 && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}


	void BCNoSlipMap::execute()
	{
		kernel->compute();
	}		

	BCVelocityMap::BCVelocityMap(SPLBGK nm, 
	                             SPPositionFunction v,
	                             SPAbstractDataWithGhostNodes map):
		BCondWithMap(map,nm->vectorTemplate),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),//< Important _BASIC has better performance
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),	
		num(nm),
		velocity(v)
	{	
	}

	BCVelocityMap::BCVelocityMap(SPLBGK nm, 
	                             SPPositionFunction v,
	                             SPAbstractDataWithGhostNodes map,
	                             SPAbstractDataWithGhostNodes computationalDomain):
		BCondWithMap(map,computationalDomain,nm->vectorTemplate),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),//< Important _BASIC has better performance
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		velocity(v)
	{	
	}
		
	BCVelocityMap::~BCVelocityMap()
	{
	}


	void BCVelocityMap::init()
	{		
		acl::ExpressionContainer kk;

		auto fX(generateDCFullSafe(num->getF()));

		initMapInfrastructure(kk);	
		acl::TypeID type(getElementType(num->getF()->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		auto isBoundary(generateVEPrivateVariable(1,typeI));

		unsigned int nDir(mapTVE->vectorTemplate->vectors.size());
		auto & block(fX->getBlock());
		unsigned int nDim(nD(block));
		asl::Index2PositionACL i2p(block,type);
		kk << i2p.initPosition;

		auto coef(generateVEPrivateVariable(1,type));
		auto rho(generateVEPrivateVariable(1,type));		
		auto mom(generateVEPrivateVariable(nDim,type));

		kk  << (coef = generateVEConstant(0.)) 
		    << (rho = generateVEConstant(0.))
			<< (mom = generateVEConstantN(nDim,0.));
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isComputationNode(i) );
			int directionShift(block.c2i(templ->vectors[i]));
			auto vBulk(subVE(fX->getEContainer(), 
			                 templ->invertVectors[i]));			                                       
			vBulk[0] = acl::generateShiftedElement(vBulk[0],directionShift);
			auto w(templ->quasiparticlesCoefs[i]);
			kk << (coef += select(generateVEConstant(0.), generateVEConstant(w), isBoundary, type));
			kk << (rho += select(generateVEConstant(0.), vBulk * w, isBoundary, type));
		}
		kk << (rho/=coef);
		kk << (mom = velocity->value(i2p.position)*rho); 
		
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isGhostNode(0) && isComputationNode(i) );

			auto vGhost(subVE(fX->getEContainer(), i));
			kk << (acl::assignmentSafe(vGhost, 
			                           select(vGhost,
			                                  rho +  mom * templ->vectors[i] * 3., 
			                                  isBoundary, 
			                                  type)));
		}

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((map->getEContainer() <= 0 && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}


	void BCVelocityMap::execute()
	{
		kernel->compute();
	}		

	BCConstantPressureVelocityMap::BCConstantPressureVelocityMap(SPLBGK nm, 
			          	           acl::VectorOfElements p, 
	                               SPAbstractDataWithGhostNodes map):
		BCondWithMap(map,nm->vectorTemplate),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),//< Important _BASIC has better performance
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		pressure(p),
		velocity(generateVEConstantN(nD(map->getBlock()),0.))
	{	
	}

	BCConstantPressureVelocityMap::BCConstantPressureVelocityMap(SPLBGK nm, 
			          	         acl::VectorOfElements p,
	                             acl::VectorOfElements v,
	                             SPAbstractDataWithGhostNodes map):
		BCondWithMap(map,nm->vectorTemplate),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),//< Important _BASIC has better performance
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		pressure(p),
		velocity(v)
	{	
	}
		
	BCConstantPressureVelocityMap::~BCConstantPressureVelocityMap()
	{
	}


	void BCConstantPressureVelocityMap::init()
	{		
		acl::ExpressionContainer kk;

		auto fX(generateDCFullSafe(num->getF()));

		initMapInfrastructure(kk);	
		acl::TypeID type(getElementType(num->getF()->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		auto isBoundary(generateVEPrivateVariable(1,typeI));

		unsigned int nDir(mapTVE->vectorTemplate->vectors.size());
		
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isGhostNode(0) && isComputationNode(i) );
			acl::VectorOfElements fEq(1);
			if (num->getCompressible())
				acl::copy(pressure  + 
				          3. * pressure * (velocity * templ->vectors[i]) +
				          4.5 * pressure * (velocity * templ->vectors[i]) * (velocity * templ->vectors[i]) -
				          1.5 * pressure * (velocity*velocity),
				          fEq);
			else
				acl::copy(pressure + 
				          3. * velocity * templ->vectors[i]+
				          4.5 * (velocity * templ->vectors[i]) * (velocity * templ->vectors[i]) -
				          1.5 * (velocity*velocity),
				          fEq);
			
			auto vGhost(subVE(fX->getEContainer(), i));
			kk << (acl::assignmentSafe(vGhost, 
			                           select(vGhost,
			                                  fEq, 
			                                  isBoundary, 
			                                  type)));
		}

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((map->getEContainer() <= 0 && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}


	void BCConstantPressureVelocityMap::execute()
	{
		kernel->compute();
	}		


	BCTransportLimitedDepositionMap::
		BCTransportLimitedDepositionMap(SPLBGK nm, 
		                                acl::VectorOfElements p,
		                                acl::VectorOfElements lF,
		                                SPAbstractDataWithGhostNodes map):
		BCondWithMap(map,nm->vectorTemplate),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),//< Important _BASIC has better performance
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),	
		num(nm),
		p0(p),
		limitingFactor(lF)
	{	
	}
		
	BCTransportLimitedDepositionMap::~BCTransportLimitedDepositionMap()
	{
	}


	void BCTransportLimitedDepositionMap::init()
	{		
		acl::ExpressionContainer kk;

		auto fX(generateDCFullSafe(num->getF()));

		initMapInfrastructure(kk);	
		acl::TypeID type(getElementType(num->getF()->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		auto isBoundary(generateVEPrivateVariable(1,typeI));

		unsigned int nDir(mapTVE->vectorTemplate->vectors.size());
		auto & block(fX->getBlock());
		asl::Index2PositionACL i2p(block,type);
		kk << i2p.initPosition;

		auto coef(generateVEPrivateVariable(1,type));
		auto rho(generateVEPrivateVariable(1,type));		

		kk  << (coef = generateVEConstant(0.)) 
		    << (rho = generateVEConstant(0.));

		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isComputationNode(i) );
			int directionShift(block.c2i(templ->vectors[i]));
			auto vBulk(subVE(fX->getEContainer(), 
			                 templ->invertVectors[i]));			                                       
			vBulk[0] = acl::generateShiftedElement(vBulk[0],directionShift);
			auto w(templ->quasiparticlesCoefs[i]);
			kk << (coef += select(generateVEConstant(0.), generateVEConstant(w), isBoundary, type));
			kk << (rho += select(generateVEConstant(0.), vBulk * w, isBoundary, type));
		}
		
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isGhostNode(0) && isComputationNode(i) );

			auto vGhost(subVE(fX->getEContainer(), i));
			auto jSubl(rho-(rho-coef*p0)/(1.+coef*limitingFactor));
			kk << (acl::assignmentSafe(vGhost, 
			                           select(vGhost,
			                                  jSubl/coef, 
			                                  isBoundary, 
			                                  type)));
		}

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((map->getEContainer() <= 0 && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}


	void BCTransportLimitedDepositionMap::execute()
	{
		kernel->compute();
	}		

	BCKineticsLimitedDepositionMap::
		BCKineticsLimitedDepositionMap(SPLBGK nm, 
		                               acl::VectorOfElements p,
		                               acl::VectorOfElements lF,
		                               acl::VectorOfElements b,
		                               SPAbstractDataWithGhostNodes map):
		BCondWithMap(map,nm->vectorTemplate),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),//< Important _BASIC has better performance
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),	
		num(nm),
		p0(p),
		limitingFactor(lF),
		beta(b)
	{	
	}
		
	BCKineticsLimitedDepositionMap::~BCKineticsLimitedDepositionMap()
	{
	}


	void BCKineticsLimitedDepositionMap::init()
	{		
		acl::ExpressionContainer kk;

		auto fX(generateDCFullSafe(num->getF()));

		initMapInfrastructure(kk);	
		acl::TypeID type(getElementType(num->getF()->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		auto isBoundary(generateVEPrivateVariable(1,typeI));

		unsigned int nDir(mapTVE->vectorTemplate->vectors.size());
		auto & block(fX->getBlock());
		asl::Index2PositionACL i2p(block,type);
		kk << i2p.initPosition;

		auto coef(generateVEPrivateVariable(1,type));
		auto rho(generateVEPrivateVariable(1,type));		

		kk  << (coef = generateVEConstant(0.)) 
		    << (rho = generateVEConstant(0.));

		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isComputationNode(i) );
			int directionShift(block.c2i(templ->vectors[i]));
			auto vBulk(subVE(fX->getEContainer(), 
			                 templ->invertVectors[i]));			                                       
			vBulk[0] = acl::generateShiftedElement(vBulk[0],directionShift);
			auto w(templ->quasiparticlesCoefs[i]);
			kk << (coef += select(generateVEConstant(0.), generateVEConstant(w), isBoundary, type));
			kk << (rho += select(generateVEConstant(0.), vBulk * w, isBoundary, type));
		}
		
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isGhostNode(0) && isComputationNode(i) );

			auto vGhost(subVE(fX->getEContainer(), i));
			auto jSubl(rho-beta*(rho/coef-p0)/(1.+limitingFactor));
			kk << (acl::assignmentSafe(vGhost, 
			                           select(vGhost,
			                                  jSubl/coef, 
			                                  isBoundary, 
			                                  type)));
		}

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((map->getEContainer() <= 0 && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}


	void BCKineticsLimitedDepositionMap::execute()
	{
		kernel->compute();
	}		

	ComputeSurfaceFluxMap::ComputeSurfaceFluxMap(SPLBGK nm,
	                                                 SPDataWithGhostNodesACLData fF,
			                                         SPAbstractDataWithGhostNodes map):
		BCondWithMap(map,nm->vectorTemplate),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)), //< Important _BASIC has better performance
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		fluxField(fF)		
	{	
	}

	ComputeSurfaceFluxMap::~ComputeSurfaceFluxMap()
	{
	}


	void ComputeSurfaceFluxMap::init()
	{		
		acl::ExpressionContainer kk;

		auto fX(generateDCFullSafe(num->getF()));

		initMapInfrastructure(kk);	
		acl::TypeID type(getElementType(num->getF()->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		auto isBoundary0(generateVEPrivateVariable(1,typeI));
		auto isBoundaryI(generateVEPrivateVariable(1,typeI));

		auto flux(generateVEPrivateVariable(1,type));
		
		unsigned int nDir(templ->vectors.size());
		auto & block(fX->getBlock());

		kk << (flux = acl::generateVEConstant(0));
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary0 = isGhostNode(0) && isComputationNode(i) );
			kk << (isBoundaryI = isGhostNode(i) && isComputationNode(0) );
			
			auto f0(subVE(fX->getEContainer(), i));
			int directionShift(block.c2i(templ->vectors[i]));
			auto fI(subVE(fX->getEContainer(), 
			              templ->invertVectors[i]));
			fI[0] = acl::generateShiftedElement(fI[0],directionShift);
			
			kk << (flux += select((f0-fI)*templ->quasiparticlesCoefs[i], isBoundary0, type));
			kk << (flux += select((fI-f0)*templ->quasiparticlesCoefs[i], isBoundaryI, type));
		}

		kk << (acl::assignmentSafe(fluxField->getEContainer(),flux));

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((fabs(map->getEContainer()) < .9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}


	void ComputeSurfaceFluxMap::execute()
	{
		kernel->compute();
	}


	ComputeSurfaceForceMap::ComputeSurfaceForceMap(SPLBGK nm,
	                                                 SPDataWithGhostNodesACLData fF,
			                                         SPAbstractDataWithGhostNodes map):
		BCondWithMap(map,nm->vectorTemplate),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),//< Important _BASIC has better performance
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		num(nm),
		forceField(fF)		
	{	
	}


	ComputeSurfaceForceMap::~ComputeSurfaceForceMap()
	{
	}


	void ComputeSurfaceForceMap::init()
	{		
		acl::ExpressionContainer kk;

		auto fX(generateDCFullSafe(num->getF()));

		initMapInfrastructure(kk);	
		acl::TypeID type(getElementType(num->getF()->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		auto isBoundary0(generateVEPrivateVariable(1,typeI));
		auto isBoundaryI(generateVEPrivateVariable(1,typeI));

		unsigned int nDir(templ->vectors.size());
		auto & block(fX->getBlock());
		unsigned int nDim(nD(block));

		
		auto force(generateVEPrivateVariable(nDim,type));
		

		kk<< (force = acl::generateVEConstantN(nDim,0));
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary0 = isGhostNode(0) && isComputationNode(i) );
			kk << (isBoundaryI = isGhostNode(i) && isComputationNode(0) );
			
			auto f0(subVE(fX->getEContainer(), i));
			int directionShift(block.c2i(templ->vectors[i]));
			auto fI(subVE(fX->getEContainer(), 
			              templ->invertVectors[i]));
			fI[0] = acl::generateShiftedElement(fI[0],directionShift);
			auto wc((templ->quasiparticlesCoefs[i]*AVec<double>(templ->vectors[i])));
			kk << (force += select(acl::generateVEConstantN(nDim,0),( f0+fI-2.)*wc, isBoundary0, type));
			kk << (force += select(acl::generateVEConstantN(nDim,0),(-f0-fI+2.)*wc, isBoundaryI, type));
		}

		kk << (acl::assignmentSafe(forceField->getEContainer(),force));

		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((fabs(map->getEContainer()) < .9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}


	void ComputeSurfaceForceMap::execute()
	{
		kernel->compute();
	}		
		

	SPBCond generateBCNoSlip(SPLBGK nm, 
	                         const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCNoSlip>(nm));
		addSlices(*bc, sl);	
		return bc;
	}

/*
	SPBCond generateBCNoSlipVel(SPLBGK nm, 
	                            const std::vector<SlicesNames> & sl)
	{
		unsigned int nd(nD(nm->getVectorTemplate()->vectors[0]));
		
		return generateBCConstantValueMiddlePoint(nm->getVelocity(),
		                                          AVec<double>(nd,0.),
		                                          sl,
		                                          nearestNeigboursVT(nd));	
	}
*/

	SPBCond generateBCConstantVelocity(SPLBGK nm,
	                                   AVec<> v,
	                                   const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCConstantVelocity>(nm, generateVEConstant(v)));
		addSlices(*bc, sl);	
		return bc;
		
	}


	SPBCond generateBCConstantPressure(SPLBGK nm,
	                                   double p,
	                                   const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCConstantPressure>(nm, generateVEConstant(p)));
		addSlices(*bc, sl);	
		return bc;
		
	}


	SPBCond generateBCConstantPressureVelocity(SPLBGK nm,
	                                   double p,
	                                   AVec<> v,
	                                   const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCConstantPressureVelocity>(nm, 
		                                                generateVEConstant(p),
		                                                generateVEConstant(v)));
		addSlices(*bc, sl);	
		return bc;		
	}


	SPNumMethod generateBCNoSlip(SPLBGK nm, 
	                             SPAbstractDataWithGhostNodes map)
	{
		return make_shared<BCNoSlipMap>(nm,map);
	}


	SPNumMethod generateBCNoSlipVel(SPLBGK nm, 
	                                SPAbstractDataWithGhostNodes map)
	{
		unsigned int nd(nD(nm->getVelocity()->getBlock()));
		return generateBCConstantValueMiddlePoint(nm->getVelocity(), 
		                                          AVec<>(nd,0.), 
		                                          map,
		                                          nm->vectorTemplate);
	}


	SPNumMethod generateBCNoSlipRho(SPLBGK nm, 
	                                SPAbstractDataWithGhostNodes map)
	{
		unsigned int nd(nD(nm->getVelocity()->getBlock()));
		return generateBCConstantGradient(nm->getRho(), 
		                                  AVec<>(nd,0.), 
		                                  map,
		                                  nm->vectorTemplate);
	}


	SPNumMethod generateBCVelocity(SPLBGK nm, 
	                               SPPositionFunction v, 
	                               SPAbstractDataWithGhostNodes map)
	{
		return make_shared<BCVelocityMap>(nm,v,map);
	}


	SPNumMethod generateBCVelocity(SPLBGK nm, 
	                               SPPositionFunction v, 
	                               SPAbstractDataWithGhostNodes map,
	                               SPAbstractDataWithGhostNodes computationalDomain)
	{
		return make_shared<BCVelocityMap>(nm,v,map,computationalDomain);
	}


	SPNumMethod generateBCVelocityVel(SPLBGK nm, 
	                                  SPPositionFunction v, 
	                                  SPAbstractDataWithGhostNodes map)
	{
		return generateBCConstantValue(nm->getVelocity(), v, map);
	}


	SPNumMethod generateBCConstantPressure(SPLBGK nm,
	                                   double p,
	                                   SPAbstractDataWithGhostNodes map)
	{
		auto bc(make_shared<BCConstantPressureVelocityMap>(nm, generateVEConstant(p), map));
		return bc;
		
	}


	SPNumMethod generateBCConstantPressureVelocity(SPLBGK nm,
	                                   double p,
	                                   AVec<> v,
	                                   SPAbstractDataWithGhostNodes map)
	{
		auto bc(make_shared<BCConstantPressureVelocityMap>(nm, 
		                                                   generateVEConstant(p),
		                                                   generateVEConstant(v),
		                                                   map));
		return bc;		
	}


	SPNumMethod generateBCTransportLimitedDeposition(SPLBGK nm, 
	                                                 double p0,
	                                                 double limitingFactor,            
	                                                 SPAbstractDataWithGhostNodes map)	
	{
		auto bc(make_shared<BCTransportLimitedDepositionMap>(nm, 
		                                                   generateVEConstant(p0),
		                                                   generateVEConstant(limitingFactor),
		                                                   map));
		return bc;		
	}


	SPNumMethod generateBCKineticsLimitedDeposition(SPLBGK nm, 
	                                                double p0,
	                                                double limitingFactor,            
	                                                double beta,
	                                                SPAbstractDataWithGhostNodes map)	
	{
		auto bc(make_shared<BCKineticsLimitedDepositionMap>(nm,  
		                                                    generateVEConstant(p0),
		                                                    generateVEConstant(limitingFactor),
		                                                    generateVEConstant(beta),
		                                                    map));
		return bc;		
	}

		
	SPNumMethod generateComputeSurfaceFlux(SPLBGK nm, 
	                                       SPDataWithGhostNodesACLData fF, 
	                                       SPAbstractDataWithGhostNodes map)
	{
		auto a(make_shared<ComputeSurfaceFluxMap>(nm, fF, map));
		return a;
	}


	SPNumMethod generateComputeSurfaceForce(SPLBGK nm, 
	                                       SPDataWithGhostNodesACLData fF, 
	                                       SPAbstractDataWithGhostNodes map)
	{
		auto a(make_shared<ComputeSurfaceForceMap>(nm, fF, map));
		return a;
	}

} // asl