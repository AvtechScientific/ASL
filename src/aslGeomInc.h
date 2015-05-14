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


#ifndef ASLGEOMINC_H
#define ASLGEOMINC_H

#include<memory>
#include<vector>
#include<acl/aclHardware.h>

namespace acl
{
	class KernelConfiguration;
}

namespace asl {

	template <typename T> class AVec;
	class Block;

	class AbstractDataWithGhostNodes; 
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;
	
	class DistanceFunction;
	typedef std::shared_ptr<DistanceFunction> SPDistanceFunction;

	
	
	///generates sphere \ingroup Geom
	/**
	 \param r radius
	 \param c center
	 */
	SPDistanceFunction generateDFSphere(double r, const AVec<double> & c);


	///generates infinite cylinder \ingroup Geom
	/**
	 \param r radius
	 \param l orientation
	 \param c center
	 */
	SPDistanceFunction generateDFCylinderInf(double r, 
	                                      const AVec<double> & l, 
	                                      const AVec<double> & c);

	///generates cylinder \ingroup Geom
	/**
	 \param r radius
	 \param l orientation and length
	 \param c center
	 */
	SPDistanceFunction generateDFCylinder(double r, 
	                                      const AVec<double> & l, 
	                                      const AVec<double> & c);
	
	///generates cone \ingroup Geom
	/**
	 \param r radius
	 \param l orientation and length, from apex to base
	 \param a apex
	 */
	SPDistanceFunction generateDFCone(double r, 
	                                  const AVec<double> & l, 
	                                  const AVec<double> & a);

	/// generates a plane \ingroup Geom
	/// \param n vector orthogonal to the plane
	/// \param p0 point on the plane
	SPDistanceFunction generateDFPlane(const AVec<double> & n,
	                                   const AVec<double> & p0);


	/// generates infinite prism with convex polygon at its base \ingroup Geom
	/**
		 \param points list of points in 3D space. 

		 The prism axis is oriented orthogonal to the plane of the triangle defined by 
		 a center and first two points. Points of the base polygon must
		 be provided in continuous manner.
	*/
	SPDistanceFunction generateDFConvexPolygonPrism(std::vector<AVec<double>> points);


	///generates pyramid with convex polygon at its base and apex \p a \ingroup Geom
	/**
		 \param points list of points in 3D space. 

		 the pyramid base plane is defined by a center and fist two points
		 the points should be ordered counter clock rotation for an observer placed in \p a
	*/
	SPDistanceFunction generateDFConvexPolygonPyramid(std::vector<AVec<double>> points, 
	                                                  AVec<double> a);
	

	/// generates map corresponding to external (ghost) part of the block \ingroup Geom
	SPDistanceFunction generateDFInBlock(const Block & b, unsigned int nG);


	/// \ingroup Geom
	/**
		If both operands are not initialized - it leads to an error. 
	    If only one operand is not initialized - it is ignored and the second is returned.
	*/
	SPDistanceFunction operator&(SPDistanceFunction a, SPDistanceFunction b);
	/// \ingroup Geom
	/**
	 If both operands are not initialized - it leads to an error. 
	 If only one operand is not initialized - it is ignored and the second is returned.
	*/
	SPDistanceFunction operator|(SPDistanceFunction a, SPDistanceFunction b);
	/// \ingroup Geom
	SPDistanceFunction operator-(SPDistanceFunction a);
	/// \ingroup Geom
	SPDistanceFunction normalize(SPDistanceFunction a, double dx);

	/// \ingroup Geom
	SPAbstractDataWithGhostNodes generateDataContainer_SP(const Block &b, 
	                                                      SPDistanceFunction df, 
	                                                      unsigned int gN,
	                                                      acl::TypeID t);
	

	/// Initialize \p d by \p f
	void initData(SPAbstractDataWithGhostNodes d, 
	              SPDistanceFunction f);
	/// Initialize \p d by \p f
	void initData(SPAbstractDataWithGhostNodes d, 
	              SPDistanceFunction f, 
	              const acl::KernelConfiguration & k);
	
} // asl

#endif // ASLDISTANCEFUNCTION
