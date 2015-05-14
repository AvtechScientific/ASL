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


#include "aslGeomInc.h"
#include <math/aslDistanceFunction.h>
#include "acl/aclMath/aclVectorOfElements.h"
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include <data/aslBlocks.h>
#include <aslGenerators.h>
#include <math/aslIndex2Position.h>
#include <acl/Kernels/aclKernelConfigurationTemplates.h>
#include <acl/Kernels/aclKernel.h>
#include <data/aslDataWithGhostNodes.h>

namespace asl
{

	SPDistanceFunction operator&(SPDistanceFunction a, SPDistanceFunction b)
	{
		if ((a.get() == 0) && (b.get() == 0))
			errorMessage("DistanceFunction::operator& - both operands are not initialized");

		if (a.get() == 0)
			return b;

		if (b.get() == 0)
			return a;
		
		return SPDistanceFunction(new DistFOperationAnd(a, b));
	}
	
	
	SPDistanceFunction operator|(SPDistanceFunction a, SPDistanceFunction b)
	{
		if ((a.get() == 0) && (b.get() == 0))
			errorMessage("DistanceFunction::operator| - both operands are not initialized");

		if (a.get() == 0)
			return b;

		if (b.get() == 0)
			return a;

		return SPDistanceFunction(new DistFOperationOr(a, b));
	}

	
	SPDistanceFunction operator-(SPDistanceFunction a)
	{
		if (a.get() == 0)
			errorMessage("DistanceFunction::operator- - operand is not initialized");

		return SPDistanceFunction(new DistFOperationInversion(a));
	}


	SPDistanceFunction normalize(SPDistanceFunction a, double dx)
	{
		if (a.get() == 0)
			errorMessage("DistanceFunction::normalize - argument is not initialized");

		return SPDistanceFunction(new DistFNormalization(a, dx));
	}
	
	SPDistanceFunction generateDFSphere(double r, const AVec<double> & c)
	{
		return SPDistanceFunction(new DistFSphere(r,c));
	}

	SPDistanceFunction generateDFCylinderInf(double r, 
	                                      const AVec<double> & l, 
	                                      const AVec<double> & c)
	{
		return SPDistanceFunction(new DistFCylinder(r, l, c));
	}

	SPDistanceFunction generateDFCylinder(double r, 
	                                      const AVec<double> & l, 
	                                      const AVec<double> & c)
	{
		auto cyl(generateDFCylinderInf(r, l, c));
		auto p1(generateDFPlane(l, c + l * .5));
		auto p2(generateDFPlane(-l, c - l * .5));
		return cyl & p1 & p2;
	}

	SPDistanceFunction generateDFCone(double r, 
	                                  const AVec<double> & l, 
	                                  const AVec<double> & a)
	{
		auto cone(SPDistanceFunction(new DistFCone(r/sqrt(l2(l)), l, a)));
		auto p1(generateDFPlane(l, a + l));
		auto p2(generateDFPlane(-l, a));
		return cone & p1 & p2;
	}
		
	SPDistanceFunction generateDFPlane(const AVec<double> & n, const AVec<double> & p0)
	{
		return SPDistanceFunction(new DistFPlane(n,p0));
	}
		
	SPDistanceFunction generateDFConvexPolygonPrism(vector<AVec<double>> points)
	{
		AVec<double> center(points[0]);
		for(unsigned int i(1); i < points.size(); ++i)
			center+=points[i];
		center = center / double(points.size());

		auto axis(crossProduct(points[0]-center, points[1]-center));

		
		auto res(generateDFPlane(crossProduct(points[1]-points[0], axis),
		                         points[0]));
		for(unsigned int i(1); i < points.size(); ++i)
			res = res & generateDFPlane(crossProduct(points[(i+1) % points.size()]-points[i], axis),
			                            points[i]);
		return res;
	}

	SPDistanceFunction generateDFConvexPolygonPyramid(vector<AVec<double>> points, 
	                                                AVec<double> a)
	{
		AVec<double> center(points[0]);
		for(unsigned int i(1); i < points.size(); ++i)
			center+=points[i];
		center = center / double(points.size());

		auto axis(crossProduct(points[0]-center, points[1]-center));

		
		auto res(generateDFPlane(crossProduct(points[1]-points[0], a - points[0]),
		                         points[0]));
		for(unsigned int i(1); i < points.size(); ++i)
			res = res & generateDFPlane(crossProduct(points[(i+1) % points.size()]-points[i], 
			                                         a - points[i]),
			                            points[i]);
		res = res & generateDFPlane(-axis,center);
		return res;
	}
		
	SPDistanceFunction generateDFInBlock(const Block & b, unsigned int nG)	
	{
		unsigned int nd(nD(b));
		AVec<> n(nd, 0.);

		AVec<> nc(n);
		nc[0] = 1.;
		AVec<> offs(nd, (double(nG) - .5) * b.dx);
		AVec<> p0(b.position + offs);
		AVec<> pE(b.getBPosition() - offs);
		auto res(generateDFPlane(nc, p0) | generateDFPlane(-nc, pE));
		if(nd >= 2 )
		{
			nc = n; nc[1] = 1;
			res = res | generateDFPlane(nc, p0) | generateDFPlane(-nc, pE);
		}
		if(nd >= 3 )
		{
			nc = n; nc[2] = 1;
			res = res | generateDFPlane(nc, p0) | generateDFPlane(-nc, pE);
		}
		
		return res;
	}

	SPAbstractDataWithGhostNodes generateDataContainer_SP(const Block &b, 
	                                                      SPDistanceFunction df, 
							      unsigned int gN,
	                                                      acl::TypeID t)
	{
		asl::Index2PositionACL i2p(offset(b, gN), t);

		return generateDataContainer_SP(b, df->getDistance(i2p.positionWithInit), gN);
	}


	void initData(SPAbstractDataWithGhostNodes d, 
	              SPDistanceFunction f)
	{		
		initData(d, f, acl::KERNEL_SIMD);
	}
		

	void initData(SPAbstractDataWithGhostNodes d, 
	              SPDistanceFunction f, 
	              const acl::KernelConfiguration & kernelConfig)
	{
		if (f.get() == 0)
			errorMessage("DistanceFunction::initData - argument is not initialized");

		auto type(getElementType(d->getEContainer()));
		asl::Index2PositionACL i2p(d->getBlock(),type);

		acl::Kernel k(kernelConfig);
		k << i2p.initPosition;
		k << acl::assignmentSafe(d->getEContainer(), f->getDistance(i2p.position));
		k.setup();
		k.compute();
	}
		
		
} // asl
