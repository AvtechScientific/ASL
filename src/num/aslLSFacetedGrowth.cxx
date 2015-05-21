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


#include "aslLSFacetedGrowth.h"
#include <math/aslTemplates.h>
#include <data/aslDataWithGhostNodes.h>
#include "acl/acl.h"
#include "acl/aclHardware.h"
#include "acl/Kernels/aclKernel.h"
#include "acl/Kernels/aclKernelConfigurationTemplates.h"
#include "acl/Kernels/aclKernelMerger.h"
#include "acl/aclGenerators.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include "acl/aclMath/aclBarycentric.h"
#include "acl/aclMath/aclMathAlg.h"
#include <algorithm>
#include <math/aslTemplates.h>
#include <math/aslTemplateVE.h>
#include <math/aslTemplatesExtras.h>
#include <math/aslTemplateVEExtras.h>
#include <math/aslDistanceFunction.h>

#include <writers/aslWriter.h>

namespace asl
{

	using acl::VectorOfElements;
	using acl::generateVEConstant;
	using acl::generateVEConstantN;
	
	inline bool isIn(int i, AVec<int> a)
	{
		bool b(false);
		for(unsigned int j(0); j < a.getSize(); ++j )
			b|=(a[j]==i);
		return b;
	}

	CrystallographicParameters::CrystallographicParameters (const vector<AVec<>> & dir, 
		                                                  const vector<double> & bSt,
		                                                  const vector<double> & bDisl,
		                                                  const double bRough):
		directions(dir),
		betaSt(bSt),
		betaDisl(bDisl),
		betaRough(bRough)
	{
	}
		
	CrystallographicParameters::CrystallographicParameters (const double bRough):
		betaRough(bRough)
	{
	}
		
	CrystallographicParameters::CrystallographicParameters (){}
	

	void CrystallographicParameters::directionCode(VectorOfElements normal,
		                                           VectorOfElements direction,
	                                               VectorOfElements cosTheta,
		                                           acl::ExpressionContainer & k)
	{
		findNearestDirectionCode(directions, normal, direction, cosTheta, k);
	}

	void CrystallographicParameters::addFacet(const AVec<> normal, double bSt, double bDisl)
	{
		directions.push_back(normal);
		betaSt.push_back(bSt);
		betaDisl.push_back(bDisl);		
	}
	

	void CrystallographicParameters::init(acl::TypeID type)
	{
		copy(acl::generateVEPrivateArray(directions, type), directionsACL);
		copy(acl::generateVEPrivateArray(betaSt, type), betaStACL);
		copy(acl::generateVEPrivateArray(betaDisl, type), betaDislACL);	
	}

	acl::VectorOfElements 
			CrystallographicParameters::velocity(acl::VectorOfElements supersaturation, 
		                                         acl::VectorOfElements dir,
	                                             acl::VectorOfElements sinTheta)
	{
		acl::TypeID type(getElementType(supersaturation));

		auto bSt(acl::excerpt(betaStACL,dir));
		auto bDisl(acl::excerpt(betaDislACL,dir));
		auto bR(generateVEConstant(betaRough));
		return supersaturation*acl::min(bR, 
		                                acl::max(sinTheta * bSt,
		                                         supersaturation*bDisl,
		                                         type),
		                                type);
	}

	acl::VectorOfElements 
			CrystallographicParameters::velocity(const acl::VectorOfElements & supersaturation, 
		                                         const acl::VectorOfElements & dir,
	                                             const acl::VectorOfElements & sinTheta,
			                                     const acl::VectorOfElements & stepVelocityLimit)
	{
		acl::TypeID type(getElementType(supersaturation));

		auto vSt(supersaturation*acl::excerpt(betaStACL,dir)*sinTheta);
		auto vDisl(supersaturation*supersaturation*acl::excerpt(betaDislACL,dir));
		auto vR(supersaturation*generateVEConstant(betaRough));
		return acl::min(vR, acl::max(min(vSt, sign(stepVelocityLimit), type), vDisl, type), type);
	}

		
	VectorOfElements 
			CrystallographicParameters::stepFactor(const VectorOfElements & dir,
			                                       const VectorOfElements & position)
	{
		auto dirVec(acl::excerpt(directionsACL,dir));
		return dirVec*position;
	}

		
	
	LSFacetedGrowth::~LSFacetedGrowth()
	{}
		
	LSFacetedGrowth::LSFacetedGrowth()
	{
	}
	
	LSFacetedGrowth::LSFacetedGrowth(Data df, DataGen c):
		LevelSetLinear(df),
		superSaturation(c)
	{		
		// this line is necessery since the kernel is initialized in a perent class by KERNEL_SIMDUA 
		kernel.reset(new acl::Kernel(acl::KERNEL_BASIC)); 
	}	
		
	void LSFacetedGrowth::initVelocityComputation()
	{
		acl::ExpressionContainer kk;
//		auto & kk(*kernel);
		acl::TypeID type(getElementType(distanceField->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);
		unsigned int nd(nD(*vectorTemplate));		
		unsigned int nv(vectorTemplate->vectors.size());
		//number of crystallographyc directions
		unsigned int ncd(crystallography.directions.size());

		crystallography.init(type);
		
		TemplateVE sSatT(*superSaturation, *vectorTemplate);
		kk << sSatT.initValues;


		auto isBoundary(generateVEPrivateVariable(1,typeI));
		auto normal(generateVEPrivateVariable(nd,type));
		auto velC(generateVEPrivateVariable(1,type));
		auto counter(generateVEPrivateVariable(1, type));
		auto sSatNode(generateVEPrivateVariable(1, type));
		
		auto vto(distanceTVE->vto);
		unsigned int nCells(vto->elementaryCells.size());

		vector<int> vectorInit(ncd,0);
		auto directionsN(generateVEPrivateArray(vectorInit,typeI));
		auto velocitiesN(generateVEPrivateArray(vectorInit,type));		
		
		vector<acl::VectorOfElements> normalC;
		kk << gcGradientAllCells(*distanceTVE,  normalC);
		vector<acl::VectorOfElements> directionC(nCells);
		vector<acl::VectorOfElements> thetaC(nCells);
		vector<acl::VectorOfElements> stepFactorC(nCells);
		for(unsigned int j(0); j < nCells; ++j)
		{
			kk << gcNormalize(normalC[j]);
			copy(generateVEPrivateVariable(1,typeI), directionC[j]);
			copy(generateVEPrivateVariable(1,type), thetaC[j]);
			crystallography.directionCode(normalC[j], directionC[j], thetaC[j], kk);
			kk << (thetaC[j] = sqrt(fabs(1.- thetaC[j]*thetaC[j])));

			copy(generateVEPrivateVariable(1,type), stepFactorC[j]);
			kk<< (stepFactorC[j] = crystallography.stepFactor(directionC[j],
			                                                  getBoundaryCenter(j)));
		}

		kk << (counter = generateVEConstant(0.));
		lVelocities.resize(nv);
		copy(generateVEPrivateVariable(nd,type), lVelocities[0]);
		kk << (lVelocities[0] = generateVEConstantN(nd,0.));
		
		for(unsigned int i(1); i < nv; ++i)
		{
			copy(generateVEPrivateVariable(nd,type), lVelocities[i]);
			kk << (normal = generateVEConstantN(nd,0.));

			kk << (sSatNode = getValueOnBoundary(sSatT.values, i));
			for(unsigned int j(0); j < nCells; ++j)
				if(isIn(i, vto->elementaryCells[j]))
				{
					kk << (velC = crystallography.velocity(sSatNode, 
					                                       directionC[j], 
					                                       thetaC[j],
					                                       stepFactorC[j] - 
					                                       crystallography.stepFactor(directionC[j], 
					                                                                  getBoundaryPoint(i))));
					auto velL(excerpt(velocitiesN,directionC[j]));
					auto dirL(excerpt(directionsN,directionC[j]));
					kk << (velL = select(velC, max(velL, velC), dirL));
					kk << (dirL = generateVEConstant(1));
				}
			kk << (normalC[0] =	generateVEConstantN(nd,0));
			for(unsigned int j(0); j < ncd; ++j)
			{
				auto ind(generateVEConstant(j));
				kk << (normalC[0] += select(generateVEConstantN(nd,0),
				                         generateVEConstant(crystallography.directions[j]),
				                         excerpt(directionsN,ind), 
				                         type)); 
			}
			kk << (thetaC[0] = max(sqrt(l2(normalC[0])), 
			                       generateVEConstant(0.1), 
			                       type));
			kk << (normalC[0] = normalC[0] / thetaC[0]);
			//.5 is a large value
			kk << (velC = generateVEConstant(.5));
			for(unsigned int j(0); j < ncd; ++j)
			{
				auto ind(generateVEConstant(j));
				auto dirVec(generateVEConstant(crystallography.directions[j]));
				auto velL(excerpt(velocitiesN,ind));
				auto dirL(excerpt(directionsN,ind));
				kk << (velC = select(velC,
				                     min (velC, 
				                          velL /*/ max((normalC[0] * dirVec ),
				                                     generateVEConstant(0.01),
				                                     type)*/, 
				                          type),
				                     dirL,
				                     type));
			}

			kk << (isBoundary = isBoundaryDir(i));
			
			kk << (lVelocities[i] = select(acl::generateVEConstantN(nd,0.),
			                               normalC[0] * velC ,
			                               isBoundary,
			                               type));
			kk << (counter += select(generateVEConstant(0.), 
			                         generateVEConstant(1.),
			                         isBoundary,
			                         type)); 
			kk << (lVelocities[0] += lVelocities[i]);
		}   
		kk << (counter = max(counter, generateVEConstant(1.), type));
		kk << (lVelocities[0] = lVelocities[0] / counter );

		kernel->addExpression(acl::elementOperators::
			                      ifElse(acl::differentSign(distanceTVE->values)[0],
		                                 kk.expression, {} ));
		Writer::current->addVector("vel", lVelocities[0], *kernel);		
	}
} // asl

