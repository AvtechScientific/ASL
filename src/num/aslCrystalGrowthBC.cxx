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


#include "aslCrystalGrowthBC.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include "acl/acl.h"
#include "acl/aclHardware.h"
#include "acl/aclGenerators.h"
#include <acl/Kernels/aclKernel.h>
#include "acl/Kernels/aclKernelConfigurationTemplates.h"
#include <aslGenerators.h>
#include <math/aslTemplateVE.h>
#include "acl/aclMath/aclMathAlg.h"
#include <math/aslIndex2Position.h>
#include <math/aslPositionFunction.h>
#include <math/aslDistanceFunctionAlg.h>

namespace asl
{

	BCLinearGrowthMap::BCLinearGrowthMap(Data d, 
	                                     const acl::VectorOfElements & ceq, 
	                                     const acl::VectorOfElements & b, 
	                                     Data map,
	                                     const VectorTemplate *const t):
		BCondWithMap(map, t),
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		data(d),
		cEq(ceq),
		beta(b)
	{
	}

	BCLinearGrowthMap::BCLinearGrowthMap(Data d, 
	                                     const acl::VectorOfElements & ceq, 
	                                     const acl::VectorOfElements & b, 
	                                     Data map,
	                                     Data computationalDomain,
	                                     const VectorTemplate *const t):
		BCondWithMap(map, computationalDomain, t),
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		data(d),
		cEq(ceq),
		beta(b)
	{
	}
		
	BCLinearGrowthMap::~BCLinearGrowthMap()
	{
	}

	void BCLinearGrowthMap::init()
	{		
		acl::ExpressionContainer kk;

		unsigned int nd(nD(data->getBlock()));
		unsigned int nc(data->getEContainer().size());
		
		auto dataX(generateDCFullSafe(data));
		
		initMapInfrastructure(kk);
					
		acl::TypeID type(getElementType(data->getEContainer()));

		vector<TemplateVE> dataTVE(nc);
		for(unsigned int ii(0); ii < nc; ++ii)
		{
			dataTVE[ii].init(*dataX, *templ, ii);
			kk << dataTVE[ii].initValues;
		}

		auto normal(generateVEPrivateVariable(nd,type));		
		auto vLocal(generateVEPrivateVariable(nc,type));
		auto counter(generateVEPrivateVariable(1,type));
		auto dirCoef(generateVEPrivateVariable(1,type));
		
		unsigned int nDir(templ->vectors.size());
		
		kk << (vLocal = acl::generateVEConstantN(nc,0.));
		kk << (counter = acl::generateVEConstant(0.));
		kk << (normal = gradient(*mapTVE));
		kk << (gcNormalize(normal));
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (dirCoef = (normal * normalize(templ->vectors[i])) + 1.);
			kk << (dirCoef *= dirCoef );
			auto xx(exBoundaryX(*mapTVE,i));
			for(unsigned int ii(0); ii < nc; ++ii)
			{
				auto vBulk(subVE(dataTVE[ii].values,i));
				auto anv(subVE(beta,ii)*(normal * templ->vectors[i]));
				auto c(subVE(cEq,ii));
				kk << (subVE(vLocal,ii) += select(acl::generateVEConstant(0.), 
					    				          (vBulk * (1.+ anv * xx) - anv*c) / 
				                                  (1. - anv*(1. - xx)), 
			    		           			      isComputationNode(i),
			        		       				  type) * dirCoef);
			}
			kk << (counter+= select(acl::generateVEConstant( 0.), 
			                        dirCoef, 
			                        isComputationNode(i),
			                        type));
		}	
		kk << (vLocal /= max(counter,acl::generateVEConstant(.01)));
		for(unsigned int ii(0); ii < nc; ++ii)
		{
			kk << (acl::assignmentSafe(subVE(data->getEContainer(),ii), 
		                               select(subVE(dataTVE[ii].values,0),
		                                      subVE(vLocal,ii),
		                                      isGhostNode(0), 
		                                      type)));
		}
		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((map->getEContainer() <= 0 && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}

	void BCLinearGrowthMap::execute()
	{
		kernel->compute();
	}

	BCLinearGrowthMap1::BCLinearGrowthMap1(Data d, 
	                                     const acl::VectorOfElements & ceq, 
	                                     const acl::VectorOfElements & b, 
	                                     Data map,
	                                     const VectorTemplate *const t):
		BCondWithMap(map, t),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		data(d),
		cEq(ceq),
		beta(b)
	{
	}

	BCLinearGrowthMap1::BCLinearGrowthMap1(Data d, 
	                                     const acl::VectorOfElements & ceq, 
	                                     const acl::VectorOfElements & b, 
	                                     Data map,
	                                     Data computationalDomain,
	                                     const VectorTemplate *const t):
		BCondWithMap(map, computationalDomain, t),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		data(d),
		cEq(ceq),
		beta(b)
	{
	}
		
	BCLinearGrowthMap1::~BCLinearGrowthMap1()
	{
	}


	void BCLinearGrowthMap1::init()
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
		kk << (sna = acl::generateVEConstant(0.));
		kk << (sna1mx = acl::generateVEConstant(0.));
		kk << (snaxc = acl::generateVEConstant(0.));
		kk << (normal = gradient(*mapTVE));
		kk << (gcNormalize(normal));
		kk << gcBoundaryAreaPerGhostPoint(*mapTVE, area); 

		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isComputationNode(i));
			kk << (sw += select(acl::generateVEConstant(templ->laplasCoefs[i]), 
			                    isBoundary,type));
			kk << (swc += select(dataTVE.getValue(i) * templ->laplasCoefs[i], 
			                    isBoundary,type));
			auto a(normalize(templ->vectors[i]));
			kk << (sna += select(normal * a, 
			                     isBoundary,type));
			kk << (sna1mx += select((normal * a)*(1.-exBoundaryX(*mapTVE,i)), 
			                        isBoundary,type));
			kk << (snaxc += select((normal * a)*exBoundaryX(*mapTVE,i)*dataTVE.getValue(i), 
			                       isBoundary,type));

		}
		auto vLocal((swc*sna-area*beta*snaxc+area*beta*cEq*sna)/(area*beta*sna1mx+sw*sna));
		kk << (acl::assignmentSafe(data->getEContainer(), 
		                           select(dataTVE.getValue(0), vLocal, isGhostNode(0), type)));
		
		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((isGhostNode() && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}

	void BCLinearGrowthMap1::execute()
	{
		kernel->compute();
	}

	BCLinearGrowthMap2::BCLinearGrowthMap2(Data d, 
	                                       const acl::VectorOfElements & ceq, 
	                                       const acl::VectorOfElements & b, 
	                                       Data map,
	                                       const VectorTemplate *const t):
		BCondWithMap(map, t),
//		kernelCN(new acl::Kernel(acl::KERNEL_SIMDUA)),
//		kernelGN(new acl::Kernel(acl::KERNEL_SIMDUA)),
		kernelCN(new acl::Kernel(acl::KERNEL_BASIC)),
		kernelGN(new acl::Kernel(acl::KERNEL_BASIC)),
		data(d),
		cEq(ceq),
		beta(b)
	{
	}

	BCLinearGrowthMap2::BCLinearGrowthMap2(Data d, 
	                                       const acl::VectorOfElements & ceq, 
	                                       const acl::VectorOfElements & b, 
	                                       Data map,
	                                       Data computationalDomain,
	                                       const VectorTemplate *const t):
		BCondWithMap(map, computationalDomain, t),
//		kernelCN(new acl::Kernel(acl::KERNEL_SIMDUA)),
//		kernelGN(new acl::Kernel(acl::KERNEL_SIMDUA)),
		kernelCN(new acl::Kernel(acl::KERNEL_BASIC)),
		kernelGN(new acl::Kernel(acl::KERNEL_BASIC)),
		data(d),
		cEq(ceq),
		beta(b)
	{
	}
		
	BCLinearGrowthMap2::~BCLinearGrowthMap2()
	{
	}

	void BCLinearGrowthMap2::init()
	{		
		unsigned int nd(nD(data->getBlock()));
		unsigned int nc(data->getEContainer().size());
		unsigned int nDir(templ->vectors.size());
		acl::TypeID type(getElementType(data->getEContainer()));
		
		auto normal(generateVEPrivateVariable(nd,type));		
		auto vLocal(generateVEPrivateVariable(nc,type));
		auto counter(generateVEPrivateVariable(1,type));
		auto dirCoef(generateVEPrivateVariable(1,type));
		auto xx(generateVEPrivateVariable(1,type));

		auto dataX(generateDCFullSafe(data));
		TemplateVE dataTVE;	


		acl::ExpressionContainer kkCN;
		initMapInfrastructure(kkCN);
		
		kkCN << (vLocal = acl::generateVEConstantN(nc,0.));
		kkCN << (counter = acl::generateVEConstant(0.));
		kkCN << (normal = gradient(*mapTVE));
		kkCN << (gcNormalize(normal));
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
				kkCN << (dirCoef *= dirCoef);
				auto anv(subVE(beta,ii)*(-normal * templ->vectors[i]));
				kkCN << (vLocal += select((dataTVE.getValue(i)*(1.+anv*xx) + anv*subVE(cEq,ii)) /
				                          (1.+ anv*(1.+xx)) * dirCoef, 
				                          isGhostNode(templ->invertVectors[i]),
				                          type));
			}
			kkCN << (acl::assignmentSafe(subVE(data->getEContainer(),ii), 
			                             select(subVE(dataTVE.values,0), 
			                                    xx*dataTVE.getValue(0) + (1.-xx)*vLocal/counter, 
			                                    isComputationNode(0), 
			                                    type)));
		}
		kernelCN->addExpression(acl::elementOperators::
			                        ifElse(acl::elementOperators::any((isGhostNode() && 
			                                                           map->getEContainer() > -.9999)[0]),
			                               kkCN.expression, 
			                               {} ));

		kernelCN->setup();
		
		acl::ExpressionContainer kkGN;
		initMapInfrastructure(kkGN);
		
		kkGN << (vLocal = acl::generateVEConstantN(nc,0.));
		kkGN << (counter = acl::generateVEConstant(0.));
		kkGN << (normal = gradient(*mapTVE));
		kkGN << (gcNormalize(normal));
		kkGN << (xx = acl::generateVEConstant(0.));
		for(unsigned int i(1); i < nDir; ++i)
		{
			kkGN << (dirCoef = (normal * normalize(templ->vectors[i])) + 1.);
			kkGN << (dirCoef *= dirCoef );			
			kkGN << (counter+= select(dirCoef, isComputationNode(i), type));
			kkGN << (xx+= select(dirCoef * exBoundaryX(*mapTVE,i), 
			                     isComputationNode(i), type));			
		}
		kkGN << (xx=select(acl::generateVEConstant(1.),xx/counter,counter>0.001,type));
		kkGN << (counter = max(counter,acl::generateVEConstant(.001), type));

		for(unsigned int ii(0); ii < nc; ++ii)
		{
			dataTVE.init(*dataX, *templ, ii,false);
			kkGN << (vLocal = acl::generateVEConstant(0.));
			for(unsigned int i(1); i < nDir; ++i)
			{
				kkGN << (dirCoef = (normal * normalize(templ->vectors[i])) + 1.);
				kkGN << (dirCoef *= dirCoef );
				auto anv(subVE(beta,ii)*(-normal * templ->vectors[i]));
				kkGN << (vLocal += select((dataTVE.getValue(i)*(1.-anv*xx) + anv*subVE(cEq,ii)) /
				                          (1.+ anv*(1.-xx)) * dirCoef, 
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

	void BCLinearGrowthMap2::execute()
	{
		kernelCN->compute();
		kernelGN->compute();
	}

		
	SPNumMethod generateBCLinearGrowth(SPAbstractDataWithGhostNodes d, 
	                                   double cEq,
	                                   double beta,
	                                   SPAbstractDataWithGhostNodes map,
	                                   const VectorTemplate *const t)
	{
		return make_shared<BCLinearGrowthMap>
			      (d, 
			       acl::generateVEConstant(cEq), 
			       acl::generateVEConstant(beta), 
			       map, 
			       t);	
	}

	SPNumMethod generateBCLinearGrowth(SPAbstractDataWithGhostNodes d, 
	                                       double cEq,
	                                       double beta,
	                                       SPAbstractDataWithGhostNodes map,
	                                       SPAbstractDataWithGhostNodes computationalDomain,
	                                       const VectorTemplate *const t)
	{
		return make_shared<BCLinearGrowthMap1>
			      (d, 
			       acl::generateVEConstant(cEq), 
			       acl::generateVEConstant(beta), 
			       map, 
			       computationalDomain, 
			       t);	
	}

	SPNumMethod generateBCLinearGrowth2(SPAbstractDataWithGhostNodes d, 
	                                    double cEq,
	                                    double beta,
	                                    SPAbstractDataWithGhostNodes map,
	                                    const VectorTemplate *const t)
	{
		return make_shared<BCLinearGrowthMap2>
			(d, 
			 acl::generateVEConstant(cEq), 
			 acl::generateVEConstant(beta), 
			 map, 
			 t);	
	}

	SPNumMethod generateBCLinearGrowth2(SPAbstractDataWithGhostNodes d, 
	                                    double cEq,
	                                    double beta,
	                                    SPAbstractDataWithGhostNodes map,
	                                    SPAbstractDataWithGhostNodes computationalDomain,
	                                    const VectorTemplate *const t)
	{
		return make_shared<BCLinearGrowthMap2>
			(d, 
			 acl::generateVEConstant(cEq), 
			 acl::generateVEConstant(beta), 
			 map, 
			 computationalDomain, 
			 t);	
	}

		
} // asl

