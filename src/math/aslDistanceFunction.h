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


#ifndef ASLDISTANCEFUNCTION_H
#define ASLDISTANCEFUNCTION_H

#include <acl/aclMath/aclVectorOfElementsDef.h>
#include <aslUtilities.h>
 
namespace acl
{
	class KernelConfiguration;
}

namespace asl {

	template <typename T> class AVec;
	class Block;

	class AbstractDataWithGhostNodes; 
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;
	
	/// Abstract class allows to add in kernel distance function from a boundary of a geometrical object
	///\ingroup LevelSet
	/**
		 The class ...
	*/
	class DistanceFunction
	{
		protected:
			DistanceFunction();
		public:
			virtual ~DistanceFunction();
			virtual acl::VectorOfElements getDistance(const acl::VectorOfElements & pos)=0;
	};

	///\ingroup LevelSet
	typedef std::shared_ptr<DistanceFunction> SPDistanceFunction;	
	
	///\ingroup LevelSet
	class DistFBinaryOperation: public DistanceFunction
	{
		protected:
			SPDistanceFunction e1;
			SPDistanceFunction e2;
		public:
			DistFBinaryOperation(SPDistanceFunction a, SPDistanceFunction b);
	};

	///\ingroup LevelSet
	class DistFUnaryOperation: public DistanceFunction
	{
		protected:
			SPDistanceFunction e1;
		public:
			DistFUnaryOperation(SPDistanceFunction a);
	};
	
	///\ingroup LevelSet
	class DistFOperationAnd: public DistFBinaryOperation
	{
		public:
			DistFOperationAnd(SPDistanceFunction a, SPDistanceFunction b);
			virtual acl::VectorOfElements getDistance(const acl::VectorOfElements & pos);
	};

	///\ingroup LevelSet
	class DistFOperationOr: public DistFBinaryOperation
	{
		public:
			DistFOperationOr(SPDistanceFunction a, SPDistanceFunction b);
			virtual acl::VectorOfElements getDistance(const acl::VectorOfElements & pos);
	};

	///\ingroup LevelSet
	class DistFOperationInversion: public DistFUnaryOperation
	{
		public:
			DistFOperationInversion(SPDistanceFunction a);
			virtual acl::VectorOfElements getDistance(const acl::VectorOfElements & pos);
	};
	
	
	///\ingroup LevelSet
	class DistFSphere: public DistanceFunction
	{
		private:
			acl::VectorOfElements radius;
			acl::VectorOfElements center;
		public:
			DistFSphere(double r, const AVec<double> & c);
			virtual acl::VectorOfElements getDistance(const acl::VectorOfElements & pos);
	};

	///\ingroup LevelSet
	/**
	*/
	class DistFCylinder: public DistanceFunction
	{
		private:
			acl::VectorOfElements radius;
			acl::VectorOfElements orientation;
			acl::VectorOfElements center;
		public:
			/**
			 \param r radius
			 \param l orientation vector no normalization needed
			 \param c center
			 */
			DistFCylinder(double r, const AVec<double> & l, const AVec<double> & c);
			virtual acl::VectorOfElements getDistance(const acl::VectorOfElements & pos);
	};

	///\ingroup LevelSet
	/**
	*/
	class DistFCone: public DistanceFunction
	{
		private:
			acl::VectorOfElements tanTheta;
			acl::VectorOfElements orientation;
			acl::VectorOfElements apex;
		public:
			/**
			 \param th theta, opening angle, 0.5*aperture
			 \param l orientation vector no normalization needed
			 \param a apex position
			 */
			DistFCone(double th, const AVec<double> & l, const AVec<double> & a);
			virtual acl::VectorOfElements getDistance(const acl::VectorOfElements & pos);
	};
	
	///\ingroup LevelSet
	/**
	*/
	class DistFPlane: public DistanceFunction
	{
		private:
			acl::VectorOfElements normal;
			acl::VectorOfElements b;
		public:
			DistFPlane(AVec<double> n, AVec<double> p0);
			virtual acl::VectorOfElements getDistance(const acl::VectorOfElements & pos);
	};
	

	/// normalize so that the values are in range [1:-1]
	///\ingroup LevelSet	
	class DistFNormalization: public DistFUnaryOperation
	{
		protected:
			acl::VectorOfElements factor;
		public:
			static const double scaleFactor;
			DistFNormalization(SPDistanceFunction a, double dx);
			virtual acl::VectorOfElements getDistance(const acl::VectorOfElements & pos);
	};

	/// transforms discrete data object into continious one, DistanceFunction
	///\ingroup LevelSet	
	class DataInterpolation: public DistanceFunction
	{
		protected:
			SPAbstractDataWithGhostNodes data;
			
		public:
			DataInterpolation(SPAbstractDataWithGhostNodes d);
			virtual acl::VectorOfElements getDistance(const acl::VectorOfElements & pos);
	};
	
} // asl

#endif // ASLDISTANCEFUNCTION
