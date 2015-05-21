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


#include "aslLevelSetLinear.h"
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
#include <math/aslTemplates.h>
#include <math/aslTemplateVE.h>
#include <math/aslTemplatesExtras.h>
#include <math/aslDistanceFunction.h>


namespace asl
{

	LevelSetLinear::~LevelSetLinear()
	{}
		
	LevelSetLinear::LevelSetLinear()
	{
	}
	
	LevelSetLinear::LevelSetLinear(Data df):
		LevelSet(df)
	{		
	}	

	/// computes gradient and distance for a cell i.
	/// val is the corresponding template with values
	/// vel is velocity templates each element contains velocity vector for corresponding node		
	/// k is kernel working with the data
	/// grad is the private variable where the corresponding value is stored
	/// d is the distance (constant of the corresponding lenear function)
	/// isBoundary set of booleans corresponding to noundary existence in the current cell
	void computeLinearFunction(unsigned int i, 
	                           TemplateVE & val, 
	                           vector<acl::VectorOfElements> & vel,
                               acl::Kernel & k,
	                           const acl::VectorOfElements & grad,
	                           const acl::VectorOfElements & d,
	                           const acl::VectorOfElements & isBoundary)
	{
		typedef acl::VectorOfElements VE;
		auto vt(val.vectorTemplate);
		auto vto(vtObject(vt));
		unsigned int nv(vto->elementaryCells[0].getSize());
		VE lval(nv);
		vector<VE> lpoints(nv);
		static acl::Barycentric bac;
		                
		for(unsigned int ii(0); ii < nv; ++ii)
		{
			unsigned int j(vto->elementaryCells[i][ii]);
			lval[ii]=val.values[j];
			copy(vel[j]+vt->vectors[j], lpoints[ii]);
		}

		bac.init(lpoints);
		k<<bac.initTInv;	
		k<<(grad = bac.gradient(lval));
		k<<(d = bac.interpolate(acl::generateVEConstant(vt->vectors[0]), lval));
		k<<(isBoundary = 
		    acl::differentSign(lval) &&
		    bac.in(-d*grad/l2(grad)));
		k<<(d = select(10.*sign(subVE(val.values,0)),
		               d*rsqrt(l2(grad)),
		               isBoundary,
		               getElementType(val.values)));
	}

	/// This function compute new value of the distance function for the zerow node of the template
	/// i is the number of vector (direction)
	/// val i is the corresponding template with values
	/// vel is the corresponding velocity vector template x y and z components
			
	void computeLinearDistance(unsigned int i, 
	                           TemplateVE & val, 
	                           vector<acl::VectorOfElements> & vel,
                               acl::Kernel & k,
	                           const acl::VectorOfElements & pd,
	                           const acl::VectorOfElements & d,
	                           const acl::VectorOfElements & isBoundary)
	{
		typedef acl::VectorOfElements VE;
		auto vt(val.vectorTemplate);
		vector<VE> lval(2);
		vector<VE> lpoints(2);
		                
		copy(subVE(val.values, 0), lval[0]);
		copy(vel[0], lpoints[0]);
		copy(subVE(val.values, i), lval[1]);
		copy(vel[i]+vt->vectors[i], lpoints[1]);

		k<<(pd = lval[0]/(lval[0]-lval[1])*(lpoints[1]-lpoints[0])+lpoints[0]);
		k<<(isBoundary = acl::differentSign(acl::cat(lval[0], lval[1])));
		k<<(d = select(10.*sign(lval[0]),
		               sqrt(l2(pd))*sign(lval[0])*sign(pd*vt->vectors[i]),
		               isBoundary,
		               getElementType(val.values)));
		
	}
		
	void LevelSetLinear::initDistancesComputation()
	{
		auto & kk(*kernel);
		acl::TypeID type(getElementType(distanceField->getEContainer()));
		unsigned int nd(nD(*vectorTemplate));
		auto vto(vtObject(vectorTemplate));

		auto newDistance(generateVEPrivateVariable(1,type));
		auto newDistanceCounter(generateVEPrivateVariable(1,type));
		auto newDistanceD(generateVEPrivateVariable(1,type));
		auto newDistanceDCounter(generateVEPrivateVariable(1,type));
		auto grad(generateVEPrivateVariable(nd,type));
		auto d(generateVEPrivateVariable(1,type));
		auto isBoundary(generateVEPrivateVariable(1,acl::TYPE_SELECT[type]));
		
		kk << (newDistance = acl::generateVEConstant(0.));//large number
		kk << (newDistanceCounter = acl::generateVEConstant(0.));//large number

/*		for(unsigned int i(0); i < vto->elementaryCells.size(); ++i)
		{
			computeLinearFunction(i, *distanceTVE, lVelocities, kk, grad, d, isBoundary);
//			kk << (newDistance = minAbs(newDistance, d));
			kk << (newDistance += select(acl::generateVEConstant(0.),
			                             d, 
			                             fabs(d) < 5.,
			                             type));
			kk << (newDistanceCounter += select(acl::generateVEConstant(0.), 
			                                    acl::generateVEConstant(1.), 
			                                    fabs(d) < 5.,
			                                    type));
		}*/
		kk << (newDistanceD = acl::generateVEConstant(0.));//large number
		kk << (newDistanceDCounter = acl::generateVEConstant(0.));//large number

		for(unsigned int i(1); i < vectorTemplate->vectors.size(); ++i)
		{
			computeLinearDistance(i, *distanceTVE, lVelocities, kk, grad, d, isBoundary);
			kk << (newDistanceD = select(acl::generateVEConstant(0.),
			                             d, 
			                             fabs(d) < 5.,
			                             type));
			kk << (newDistanceDCounter += select(acl::generateVEConstant(0.), 
			                                    acl::generateVEConstant(1.), 
			                                    fabs(d) < 5.,
			                                    type));
		}
		kk << (newDistance = select(newDistanceD/max(newDistanceDCounter,
		                                             acl::generateVEConstant(.1),
		                                             type),
		                            newDistance/newDistanceCounter,
		                            newDistanceCounter > .1,
		                            type));// / 
		                     //DistFNormalization::scaleFactor);
		kk << (newDistance = select(newDistance/DistFNormalization::scaleFactor,
		                          newDistance/(1. - newDistance),
		                          newDistance*subVE(distanceTVE->values, 0) < 0));
		
//		kk << assignmentSafe (distanceFieldInternalData->getSubContainer(),  
//	               		                         select( newDistance, 
//			                                               sign(newDistance), 
//			                                               fabs(newDistance) > 1., 
//			                                               type));
		kk << assignmentSafe (distanceFieldInternalData->getSubContainer(),  
	               		                         select( newDistance, 
			                                             sign(subVE(distanceTVE->values,0)), 
			                                             newDistanceCounter + newDistanceDCounter < .1, 
			                                             type));
		
	}
				
} // asl

