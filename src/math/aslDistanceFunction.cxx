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


#include "aslDistanceFunction.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include <acl/aclGenerators.h>
#include <acl/aclHardware.h>
#include <acl/acl.h>
#include <data/aslBlocks.h>
#include <aslGenerators.h>
#include <math/aslIndex2Position.h>
#include <acl/Kernels/aclKernelConfigurationTemplates.h>
#include <acl/Kernels/aclKernel.h>
#include <data/aslDataWithGhostNodes.h>
#include <math/aslTemplateVE.h>


namespace asl
{
	DistanceFunction::DistanceFunction()
	{}	

	DistanceFunction::~DistanceFunction()
	{}	

	DistFBinaryOperation::
		DistFBinaryOperation(SPDistanceFunction a, SPDistanceFunction b):
			e1(a), e2(b)
	{
	}

	DistFUnaryOperation::
		DistFUnaryOperation(SPDistanceFunction a):
			e1(a)
	{
	}
	
	DistFOperationAnd::
		DistFOperationAnd(SPDistanceFunction a, SPDistanceFunction b):
			DistFBinaryOperation(a,b)
	{
	}
	
	acl::VectorOfElements DistFOperationAnd::getDistance(const acl::VectorOfElements & pos)
	{
		return max(e1->getDistance(pos), e2->getDistance(pos));
	}

	DistFOperationOr::
		DistFOperationOr(SPDistanceFunction a, SPDistanceFunction b):
			DistFBinaryOperation(a,b)
	{
	}
	
	acl::VectorOfElements DistFOperationOr::getDistance(const acl::VectorOfElements & pos)
	{
		return min(e1->getDistance(pos), e2->getDistance(pos));
	}

	DistFOperationInversion::
		DistFOperationInversion(SPDistanceFunction a):
			DistFUnaryOperation(a)
	{
	}

	acl::VectorOfElements DistFOperationInversion::getDistance(const acl::VectorOfElements & pos)
	{
		return -(e1->getDistance(pos));
	}
		
		
	DistFSphere::
		DistFSphere(double r, const AVec<double> & c):
			radius(acl::generateVEConstant(r)), center(acl::generateVEConstant(c))
	{
	}

	acl::VectorOfElements DistFSphere::getDistance(const acl::VectorOfElements & pos)
	{
		return sqrt(l2(pos-center))-radius;
	}
		
	DistFCylinder::
		DistFCylinder(double r, const AVec<double> & l, const AVec<double> & c):
			radius(acl::generateVEConstant(r)),
			orientation(acl::generateVEConstant(normalize(l))), 
			center(acl::generateVEConstant(c))
	{
	}

	acl::VectorOfElements DistFCylinder::getDistance(const acl::VectorOfElements & pos)
	{
		return sqrt(l2((pos-center) - ((pos-center)*orientation)*orientation))-radius;
	}

	DistFCone::
		DistFCone(double th, const AVec<double> & l, const AVec<double> & a):
			tanTheta(acl::generateVEConstant(tan(th))),
			orientation(acl::generateVEConstant(normalize(l))), 
			apex(acl::generateVEConstant(a))
	{
	}

	acl::VectorOfElements DistFCone::getDistance(const acl::VectorOfElements & pos)
	{
		return sqrt(l2((pos-apex) - ((pos-apex)*orientation)*orientation)) - 
			tanTheta*((pos-apex)*orientation);
	}
		
	acl::VectorOfElements DistFPlane::getDistance(const acl::VectorOfElements & pos)
	{
		return pos*normal-b;
	}
		
	DistFPlane::
		DistFPlane(AVec<double> n, AVec<double> p0):
			normal(acl::generateVEConstant(normalize(n))), 
		    b(acl::generateVEConstant(p0*n/sqrt(l2(n))))
	{
	}
		
	const double DistFNormalization::scaleFactor(1.8);
		
	DistFNormalization::
		DistFNormalization(SPDistanceFunction a, double dx):
			DistFUnaryOperation(a),
			factor(acl::generateVEConstant(1./dx/scaleFactor))
	{
	}
	
	acl::VectorOfElements DistFNormalization::getDistance(const acl::VectorOfElements & pos)
	{
		auto type(getElementType(pos));
		acl::VectorOfElements dist(1); 
		copy(factor*e1->getDistance(pos),dist);
		
		acl::VectorOfElements topCut(1);
		copy(acl::select(dist, acl::generateVEConstant(1.), dist>1., type),
		     topCut);
		acl::VectorOfElements bottomCut(1);
		copy(acl::select(topCut, acl::generateVEConstant(-1.), dist<-1., type),
		     bottomCut);
		return bottomCut;
	}

	acl::VectorOfElements DataInterpolation::getDistance(const acl::VectorOfElements & pos)
	{
		acl::TypeID type(getElementType(data->getEContainer()));
		auto & bl(data->getBlock());
		auto e((pos - bl.position)/bl.dx);
		auto ind(convert(acl::TYPE_SELECT[type],floor(e) * bl.c2iTransformVector,false));
		copy(e - acl::floor(e),e);

		TemplateVE dTVE;
		auto d(generateDataContainer_SP(data->getBlock(), data->getEContainer(),0));
		unsigned int nC(data->getEContainer().size());
		unsigned int nd(nD(data->getBlock()));
		acl::VectorOfElements res(nC);
		for(unsigned int i(0); i < nC; ++i)
		{
			dTVE.init(*d, *elementaryCellVT(nd), i, false);
			copy(excerpt(dTVE.values,ind),dTVE.values);
			res[i] = (interpolate(dTVE, e))[0];
		}
		return res;
	}
		
	DataInterpolation::
		DataInterpolation(SPAbstractDataWithGhostNodes d):
			data(d)
	{
	}

		
} // asl
