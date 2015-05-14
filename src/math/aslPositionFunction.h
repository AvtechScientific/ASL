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


#ifndef ASLPOSITIONFUNCTION_H
#define ASLPOSITIONFUNCTION_H

#include <acl/aclMath/aclVectorOfElementsDef.h>
#include "utilities/aslUValue.h"

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
	///\ingroup PF
	/**
		 The class ...
	*/
	class PositionFunction
	{
		protected:
			PositionFunction();
		public:
			virtual ~PositionFunction();
			virtual acl::VectorOfElements value(acl::VectorOfElements & pos)=0;
	};

	///\ingroup PF
	typedef std::shared_ptr<PositionFunction> SPPositionFunction;

	///\ingroup PF
	class PFBinaryOperation: public PositionFunction
	{
		protected:
			SPPositionFunction e1;
			SPPositionFunction e2;
		public:
			PFBinaryOperation(SPPositionFunction a, SPPositionFunction b);
	};

	///\ingroup PF
	class PFUnaryOperation: public PositionFunction
	{
		protected:
			SPPositionFunction e1;
		public:
			PFUnaryOperation(SPPositionFunction a);
	};
	
	///\ingroup PF
	class PFOperationPlus: public PFBinaryOperation
	{
		public:
			PFOperationPlus(SPPositionFunction a, SPPositionFunction b);
			virtual acl::VectorOfElements value(acl::VectorOfElements & pos);
	};

	///\ingroup PF
	class PFOperationMinus: public PFBinaryOperation
	{
		public:
			PFOperationMinus(SPPositionFunction a, SPPositionFunction b);
			virtual acl::VectorOfElements value(acl::VectorOfElements & pos);
	};

	///\ingroup PF
	class PFOperationProduct: public PFBinaryOperation
	{
		public:
			PFOperationProduct(SPPositionFunction a, SPPositionFunction b);
			virtual acl::VectorOfElements value(acl::VectorOfElements & pos);
	};
	
	///\ingroup PF
	class PFOperationInversion: public PFUnaryOperation
	{
		public:
			PFOperationInversion(SPPositionFunction a);
			virtual acl::VectorOfElements value(acl::VectorOfElements & pos);
	};
	
	///\ingroup PF
	class PFConstant: public PositionFunction
	{
		private:
			acl::VectorOfElements val;
		public:
			PFConstant(acl::VectorOfElements v);
			virtual acl::VectorOfElements value(acl::VectorOfElements & pos);
	};

	///\ingroup PF
	class PFLinear: public PositionFunction
	{
		private:
			acl::VectorOfElements gradient;
			acl::VectorOfElements b;
		public:
			PFLinear(acl::VectorOfElements g, acl::VectorOfElements p0);
			virtual acl::VectorOfElements value(acl::VectorOfElements & pos);
	};
	
	/// creates function corresponding to a velocity field \ingroup PF 
	class PFRotationField: public PositionFunction
	{
		private:
			acl::VectorOfElements axis;
			acl::VectorOfElements c;
		public:
			/**
			 \param rotationAxis length of the vector corresponds to an angular velocity value
			 */ 			
			PFRotationField(acl::VectorOfElements rotationAxis, 
			                acl::VectorOfElements center);
			virtual acl::VectorOfElements value(acl::VectorOfElements & pos);
	};

	///\ingroup PF
	class PFSign: public PFUnaryOperation
	{
		public:
			PFSign(SPPositionFunction a);
			virtual acl::VectorOfElements value(acl::VectorOfElements & pos);
	};

	///\ingroup PF
	SPPositionFunction generatePFLinear(const AVec<double> & g, double p0);
	///\ingroup PF
	template <typename T> 		
		SPPositionFunction generatePFLinear(const AVec<double> & g, UValue<T> p0);

	///\ingroup PF
	SPPositionFunction generatePFConstant(const AVec<double> & a);
	///\ingroup PF
	SPPositionFunction generatePFConstant(double a);
	///\ingroup PF
	SPPositionFunction generatePFRotationField(const AVec<double> & axis, 
	                                           const AVec<double> & c);


	///\ingroup PF
	inline SPPositionFunction operator+(SPPositionFunction a, SPPositionFunction b);
	///\ingroup PF
	inline SPPositionFunction operator-(SPPositionFunction a, SPPositionFunction b);
	///\ingroup PF
	inline SPPositionFunction operator-(SPPositionFunction a);
	///\ingroup PF
	inline SPPositionFunction operator*(SPPositionFunction a, SPPositionFunction b);
	///\ingroup PF
	inline SPPositionFunction sign(SPPositionFunction a);

	///\ingroup PF
	SPAbstractDataWithGhostNodes generateDataContainer_SP(const Block &b, 
	                                                      SPPositionFunction df, 
	                                                      unsigned int gN,
	                                                      acl::TypeID t);


	/// Initialize \p d by \p f \ingroup PF
	void initData(SPAbstractDataWithGhostNodes d, 
	              SPPositionFunction f);
	/// Initialize \p d by \p f \ingroup PF
	void initData(SPAbstractDataWithGhostNodes d, 
	              SPPositionFunction f, 
	              const acl::KernelConfiguration & k);


//---------------------------- Implementation --------------------------------


	inline SPPositionFunction operator+(SPPositionFunction a, SPPositionFunction b)
	{
		return SPPositionFunction(new PFOperationPlus(a,b));
	}
	
	inline SPPositionFunction operator-(SPPositionFunction a, SPPositionFunction b)
	{
		return SPPositionFunction(new PFOperationMinus(a,b));
	}

	inline SPPositionFunction operator-(SPPositionFunction a)
	{
		return SPPositionFunction(new PFOperationInversion(a));
	}

	inline SPPositionFunction operator*(SPPositionFunction a, SPPositionFunction b)
	{
		return SPPositionFunction(new PFOperationProduct(a,b));
	}

	inline SPPositionFunction sign(SPPositionFunction a)
	{
		return SPPositionFunction(new PFSign(a));
	}
	
}// asl

#endif // ASLPositionFunction
