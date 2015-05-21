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


#include "aslBasicBC.h"
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
#include <utilities/aslUValue.h>
#include "aslBasicBC2.h"


namespace asl
{
	BCConstantValue::BCConstantValue(Data d, const acl::VectorOfElements & v):
		BCond(d->getBlock()),
		kernel(new acl::Kernel()),
		data(d),
		value(v)
	{
	}


	void BCConstantValue::init()
	{
		loadIndicesToACL();
		kernel->clear();
		(*kernel)<<(acl::excerpt(data->getEContainer(), *indicesACL) = value);
		kernel->setup();
	}


	void BCConstantValue::execute()
	{
		kernel->compute();
	}


	void BCConstantValue::setValue(const acl::VectorOfElements & v)
	{
		 acl::copy(v,value);
	}

	BCConstantValueMap::BCConstantValueMap(Data d, 
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

		
	BCConstantValueMiddlePointMap::
		BCConstantValueMiddlePointMap(Data d, 
	                                  const acl::VectorOfElements & v, 
	                                  Data map,
	                                  const VectorTemplate *const t):
		BCondWithMap(map, t),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		data(d),
		value(v)
	{
	}

	BCConstantValueMiddlePointMap::~BCConstantValueMiddlePointMap()
	{
	}


	void BCConstantValueMiddlePointMap::init()
	{		
		acl::ExpressionContainer kk;

		unsigned int nc(data->getEContainer().size());
		
		auto dataX(generateDCFullSafe(data));

		acl::TypeID type(getElementType(data->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		initMapInfrastructure(kk);

		vector<TemplateVE> dataTVE(nc);
		for(unsigned int ii(0); ii < nc; ++ii)
		{
			dataTVE[ii].init(*dataX, *templ, ii);
			kk << dataTVE[ii].initValues;
		}

		auto vLocal(generateVEPrivateVariable(nc,type));
		auto isBoundary(generateVEPrivateVariable(1,typeI));
		auto counter(generateVEPrivateVariable(1,type));

		unsigned int nDir(mapTVE->vectorTemplate->vectors.size());
		
		kk << (vLocal = acl::generateVEConstantN(nc,0.));
		kk << (counter = acl::generateVEConstant(0.));

		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isComputationNode(i));

			for(unsigned int ii(0); ii < nc; ++ii)
			{
				auto vBulk(subVE(dataTVE[ii].values,i));			                                       
				kk << (subVE(vLocal,ii) += select(2.*subVE(value, ii) - vBulk, 
			    		           			      isBoundary,
			        		       				  type) * templ->laplasCoefs[i]);
			}
			kk << (counter+= select(acl::generateVEConstant( templ->laplasCoefs[i] ), 
			                        isBoundary,
			                        type));
		}	

		kk << (vLocal /= max(counter,acl::generateVEConstant(.01),type));
		for(unsigned int ii(0); ii < nc; ++ii)
		{
			kk << (acl::assignmentSafe(subVE(data->getEContainer(),ii), 
		                               select(subVE(dataTVE[ii].values,0), 
		                                      subVE(vLocal,ii),
		                                      isGhostNode(0), 
		                                      type)));
		}
		kernel->addExpression(acl::elementOperators::
			                      ifElse(acl::elementOperators::any((isGhostNode() && 
			                                                         map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}

	void BCConstantValueMiddlePointMap::execute()
	{
		kernel->compute();
	}

		
	BCConstantGradient::BCConstantGradient(Data d, 
		                           const acl::VectorOfElements & v, 
	                               const VectorTemplate *const t):
		BCond(d->getBlock(), t),
		kernel(new acl::Kernel()),
		data(d),
		value(v * block.dx)
	{
	}

	void BCConstantGradient::init()
	{
		loadIndicesToACL();
		loadNeighbourIndicesToACL();
		kernel->clear();
		(*kernel) << (acl::excerpt(data->getEContainer(), *indicesACL) =
		              acl::excerpt(data->getEContainer(), *neighbourIndicesACL) -
		              value);
		kernel->setup();
	}

	void BCConstantGradient::execute()
	{
		kernel->compute();
	}

	void BCConstantGradient::setValue(const acl::VectorOfElements & v)
	{
		acl::copy(v * block.dx, value);
	}

	BCConstantGradientMap::BCConstantGradientMap(Data d, 
	                                     const acl::VectorOfElements & v, 
	                                     Data map,
	                                     const VectorTemplate *const t):
		BCondWithMap(map, t),
//		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		kernel(new acl::Kernel(acl::KERNEL_BASIC)),
		data(d),
		value(v)
	{
	}

	BCConstantGradientMap::BCConstantGradientMap(Data d, 
	                                     const acl::VectorOfElements & v, 
	                                     Data map,
	                                     Data computationalDomain,
	                                     const VectorTemplate *const t):
		BCondWithMap(map, computationalDomain, t),
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		data(d),
		value(v)
	{
	}
		
	BCConstantGradientMap::~BCConstantGradientMap()
	{
	}


	void BCConstantGradientMap::init()
	{		
		acl::ExpressionContainer kk;

		unsigned int nc(data->getEContainer().size());
		unsigned int nDir(templ->vectors.size());
		unsigned int nd(nD(data->getBlock()));
		
		auto dataX(generateDCFullSafe(data));

		initMapInfrastructure(kk);
					
		acl::TypeID type(getElementType(data->getEContainer()));

		auto normal(generateVEPrivateVariable(nd,type));		
		auto vLocal(generateVEPrivateVariable(1,type));
		auto counter(generateVEPrivateVariable(1,type));
		auto dirCoef(generateVEPrivateVariable(1,type));

		kk << (counter = acl::generateVEConstant(0.));
		kk << (normal = gradient(*mapTVE));
		kk << (gcNormalize(normal));
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (dirCoef = (normal * normalize(templ->vectors[i])) + 1.);
			kk << (dirCoef *= dirCoef );			
			kk << (counter+= select(dirCoef, isComputationNode(i), type));
		}
		kk << (counter = max(counter,acl::generateVEConstant(.001), type));
		
		
		TemplateVE dataTVE;	
		for(unsigned int ii(0); ii < nc; ++ii)
		{
			dataTVE.init(*dataX, *templ, ii,false);
			kk << (vLocal = acl::generateVEConstant(0.));
			for(unsigned int i(1); i < nDir; ++i)
			{
				kk << (dirCoef = (normal * normalize(templ->vectors[i])) + 1.);
				kk << (dirCoef *= dirCoef );							
				kk << (vLocal += select((dataTVE.getValue(i) - 
				                         (normal*normalize(templ->vectors[i]))*subVE(value, ii)) *
				                        dirCoef, 
				                        isComputationNode(i),
				                        type));
			}		
			kk << (acl::assignmentSafe(subVE(data->getEContainer(),ii), 
				                       select(subVE(dataTVE.values,0), 
				                              vLocal/counter, 
				                              isGhostNode(0) && counter>0.1, 
				                              type)));
		}
		kernel->addExpression(acl::elementOperators::
		                      ifElse(acl::elementOperators::any((map->getEContainer() <= 0 && 
		                             map->getEContainer() > -.9999)[0]),
		                             kk.expression, 
		                             {} ));

		kernel->setup();
	}

/*	void BCConstantGradientMap::init()
	{		
		acl::ExpressionContainer kk;

		unsigned int nc(data->getEContainer().size());
		unsigned int nd(nD(data->getBlock()));
		
		auto dataX(generateDCFullSafe(data));

		initMapInfrastructure(kk);
					
		acl::TypeID type(getElementType(data->getEContainer()));
		acl::TypeID typeI(acl::TYPE_SELECT[type]);

		vector<TemplateVE> dataTVE(nc);
		for(unsigned int ii(0); ii < nc; ++ii)
		{
			dataTVE[ii].init(*dataX, *templ, ii);
			kk << dataTVE[ii].initValues;
		}

		auto normal(generateVEPrivateVariable(nd,type));		
		auto vLocal(generateVEPrivateVariable(nc,type));
		auto isBoundary(generateVEPrivateVariable(1,typeI));
		auto counter(generateVEPrivateVariable(1,type));
		auto dirCoef(generateVEPrivateVariable(1,type));
		
		unsigned int nDir(templ->vectors.size());
		
		kk << (vLocal = acl::generateVEConstantN(nc,0.));
		kk << (counter = acl::generateVEConstant(0.));
		kk << (normal = gradient(*mapTVE));
		kk << (gcNormalize(normal));
		for(unsigned int i(1); i < nDir; ++i)
		{
			kk << (isBoundary = isComputationNode(i));
			kk << (dirCoef = (normal * normalize(templ->vectors[i])) + 1.);
			kk << (dirCoef *= dirCoef );
			for(unsigned int ii(0); ii < nc; ++ii)
			{
				auto vBulk(subVE(dataTVE[ii].values,i));			                                       
				kk << (subVE(vLocal,ii) += select(acl::generateVEConstant(0.), 
					    				          vBulk - subVE(value, ii) * 
				                                  (normal * templ->vectors[i]), 
			    		           			      isBoundary,
			        		       				  type) * dirCoef);
			}
			kk << (counter+= select(acl::generateVEConstant( 0.), 
			                        dirCoef, 
			                        isBoundary,
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
*/
	void BCConstantGradientMap::execute()
	{
		kernel->compute();
	}

		
	BCConstantSource::BCConstantSource(Data d, cl_double v):
		BCond(d->getBlock()),
		kernel(new acl::Kernel()),
		data(d),
		value(v)
	{
	}
		
	void BCConstantSource::init()
	{
		loadIndicesToACL();
		kernel->clear();
		(*kernel)<<(acl::excerpt(data->getEContainer(),*indicesACL)+=
		         acl::generateVEVariableR(value));
		kernel->setup();
	}


	void BCConstantSource::execute()
	{
		kernel->compute();
	}


	void BCConstantSource::setValue(cl_double v)
	{
		value = v;
	}
		

	BCDirectCopier::BCDirectCopier(Data dSource, Data dDestination):
		BCondConnector(dSource->getBlock(), dDestination->getBlock()),
		kernel(new acl::Kernel()),
		source(dSource),
		destination(dDestination)
	{
	}


	void BCDirectCopier::init()
	{
		loadIndicesToACL();
		kernel->clear();
		(*kernel)<<(acl::excerpt(destination->getEContainer(),*indices2ACL)=
		         acl::excerpt(source->getEContainer(),*indices1ACL));
//		(*kernel)<<(*indices2ACL=*indices1ACL);
		kernel->setup();
	}


	void BCDirectCopier::execute()
	{
		kernel->compute();
	}

	SPBCond generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                double v, 
	                                const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCConstantValue>(d, acl::generateVEConstant(v)));
		addSlices(*bc,sl);	
		return bc;
	}		

	SPBCond generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                UValue<double> & v, 
	                                const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCConstantValue>(d, acl::generateVEVariableSP(v.p)));
		addSlices(*bc,sl);	
		return bc;		
	}
		
	SPBCond generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                UValue<AVec<float>> & v, 
	                                const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCConstantValue>(d, acl::generateVEVariableSP(v.p)));
		addSlices(*bc,sl);	
		return bc;		
	}
		
	SPBCond generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                AVec<> v, 
	                                const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCConstantValue>(d, acl::generateVEConstant(v)));
		addSlices(*bc,sl);	
		return bc;
	}
		
	SPNumMethod generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                double v, 
	                                SPAbstractDataWithGhostNodes map)
	{
		return make_shared<BCConstantValueMap>(d, acl::generateVEConstant(v), map);		
	}

	SPNumMethod generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                AVec<> v, 
	                                SPAbstractDataWithGhostNodes map)
	{
		return make_shared<BCConstantValueMap>(d, acl::generateVEConstant(v), map);		
	}

	SPNumMethod generateBCConstantValue(SPAbstractDataWithGhostNodes d, 
	                                    SPPositionFunction v, 
	                                    SPAbstractDataWithGhostNodes map)
	{
		return make_shared<BCValuePFMap>(d, v, map);		
	}
		
	SPNumMethod generateBCConstantValueMiddlePoint(SPAbstractDataWithGhostNodes d, 
	                                               double v, 
	                                               SPAbstractDataWithGhostNodes map,
	                                               const VectorTemplate *const t)
	{
		return make_shared<BCConstantValueMiddlePointMap>
			      (d, acl::generateVEConstant(v), map, t);		
	}

	SPNumMethod generateBCConstantValueMiddlePoint(SPAbstractDataWithGhostNodes d, 
	                                               AVec<> v, 
	                                               SPAbstractDataWithGhostNodes map,
	                                               const VectorTemplate *const t)
	{
		return make_shared<BCConstantValueMiddlePointMap>
			      (d, acl::generateVEConstant(v), map, t);		
	}

	SPBCond generateBCConstantGradient(SPAbstractDataWithGhostNodes d, 
	                               double v,
	                               const VectorTemplate *const t,
	                               const std::vector<SlicesNames> & sl)
	{
		auto bc(make_shared<BCConstantGradient>(d, acl::generateVEConstant(v),t));
		addSlices(*bc,sl);	
		return bc;
	}

	SPNumMethod generateBCConstantGradient(SPAbstractDataWithGhostNodes d, 
	                                       double v, 
	                                       SPAbstractDataWithGhostNodes map,
	                                       const VectorTemplate *const t)
	{
		return make_shared<BCConstantGradientMap>
			      (d, acl::generateVEConstant(v), map, t);	
	}

	SPNumMethod generateBCConstantGradient(SPAbstractDataWithGhostNodes d, 
	                                       AVec<> v, 
	                                       SPAbstractDataWithGhostNodes map,
	                                       const VectorTemplate *const t)
	{
		return make_shared<BCConstantGradientMap>
			      (d, acl::generateVEConstant(v), map, t);	
	}

	SPNumMethod generateBCConstantGradient(SPAbstractDataWithGhostNodes d, 
	                                       double v, 
	                                       SPAbstractDataWithGhostNodes map,
	                                       SPAbstractDataWithGhostNodes computationalDomain,
	                                       const VectorTemplate *const t)
	{
		return make_shared<BCConstantGradientMap>
			      (d, acl::generateVEConstant(v), map, computationalDomain, t);	
	}

		
	SPNumMethod generateBCConstantGradient2(SPAbstractDataWithGhostNodes d, 
	                                        double v, 
	                                        SPAbstractDataWithGhostNodes map,
	                                        const VectorTemplate *const t)
	{
		return make_shared<BCConstantGradientMap2>
			      (d, acl::generateVEConstant(v), map, t);	
	}

	SPNumMethod generateBCConstantGradient2(SPAbstractDataWithGhostNodes d, 
	                                        AVec<> v, 
	                                        SPAbstractDataWithGhostNodes map,
	                                        const VectorTemplate *const t)
	{
		return make_shared<BCConstantGradientMap2>
			(d, acl::generateVEConstant(v), map, t);	
	}

	SPNumMethod generateBCConstantGradient2(SPAbstractDataWithGhostNodes d, 
	                                        double v, 
	                                        SPAbstractDataWithGhostNodes map,
	                                        SPAbstractDataWithGhostNodes computationalDomain,
	                                        const VectorTemplate *const t)
	{
		return make_shared<BCConstantGradientMap2>
			(d, acl::generateVEConstant(v), map, computationalDomain, t);	
	}
		
} // asl

