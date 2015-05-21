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


#include "aslPositionFunction.h"
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
	PositionFunction::PositionFunction()
	{}	

	PositionFunction::~PositionFunction()
	{}	

	PFBinaryOperation::
		PFBinaryOperation(SPPositionFunction a, SPPositionFunction b):
			e1(a), e2(b)
	{
	}

	PFUnaryOperation::
		PFUnaryOperation(SPPositionFunction a):
			e1(a)
	{
	}
	
	PFOperationPlus::
		PFOperationPlus(SPPositionFunction a, SPPositionFunction b):
			PFBinaryOperation(a,b)
	{
	}
	
	acl::VectorOfElements PFOperationPlus::value(acl::VectorOfElements & pos)
	{
		return e1->value(pos) + e2->value(pos);
	}

	PFOperationMinus::
		PFOperationMinus(SPPositionFunction a, SPPositionFunction b):
			PFBinaryOperation(a,b)
	{
	}
	
	acl::VectorOfElements PFOperationMinus::value(acl::VectorOfElements & pos)
	{
		return e1->value(pos) - e2->value(pos);
	}

	PFOperationInversion::
		PFOperationInversion(SPPositionFunction a):
			PFUnaryOperation(a)
	{
	}

	acl::VectorOfElements PFOperationInversion::value(acl::VectorOfElements & pos)
	{
		return -(e1->value(pos));
	}
						
	PFOperationProduct::
		PFOperationProduct(SPPositionFunction a, SPPositionFunction b):
			PFBinaryOperation(a,b)
	{
	}
	
	acl::VectorOfElements PFOperationProduct::value(acl::VectorOfElements & pos)
	{
		return e1->value(pos) * e2->value(pos);
	}

	acl::VectorOfElements PFLinear::value(acl::VectorOfElements & pos)
	{
		return pos * gradient + b;
	}
		
	PFLinear::
		PFLinear(acl::VectorOfElements g, acl::VectorOfElements a):
			gradient(g), 
		    b(a)
	{
	}			
	acl::VectorOfElements PFConstant::value(acl::VectorOfElements & pos)
	{
		return val;
	}
		
	PFConstant::
		PFConstant(acl::VectorOfElements v):
			val(v)
	{
	}

	PFSign::
		PFSign(SPPositionFunction a):
			PFUnaryOperation(a)
	{
	}

	acl::VectorOfElements PFSign::value(acl::VectorOfElements & pos)
	{
		return sign(e1->value(pos));
	}
		
	SPPositionFunction generatePFLinear(const AVec<double> & g, double  p0)
	{
		return SPPositionFunction(new PFLinear(acl::generateVEConstant(g), acl::generateVEConstant(p0)));
	}

	template <typename T> 		
		SPPositionFunction generatePFLinear(const AVec<double> & g, UValue<T>  p0)
	{
		return SPPositionFunction(new PFLinear(acl::generateVEConstant(g), 
		                                       acl::generateVEVariableSP(p0.p)));
	}
	template SPPositionFunction generatePFLinear(const AVec<double> & g, UValue<double>  p0);
	template SPPositionFunction generatePFLinear(const AVec<double> & g, UValue<float>  p0);

		
	SPPositionFunction generatePFConstant(double v)
	{
		return SPPositionFunction(new PFConstant(acl::generateVEConstant(v)));
	}
		
	SPPositionFunction generatePFConstant(const AVec<double> & v)
	{
		return SPPositionFunction(new PFConstant(acl::generateVEConstant(v)));
	}
				
	SPAbstractDataWithGhostNodes generateDataContainer_SP(const Block &b, 
	                                                      SPPositionFunction df, 
	                                             		  unsigned int gN,
	                                                      acl::TypeID t)
	{
		asl::Index2PositionACL i2p(offset(b, gN), t);

		return generateDataContainer_SP(b, df->value(i2p.positionWithInit), gN);
	}

	void initData(SPAbstractDataWithGhostNodes d, 
	              SPPositionFunction f)
	{		
		initData(d, f, acl::KERNEL_SIMD);
    }
		
	void initData(SPAbstractDataWithGhostNodes d, 
	              SPPositionFunction f, 
	              const acl::KernelConfiguration & kernelConfig)
	{
		auto type(getElementType(d->getEContainer()));
		asl::Index2PositionACL i2p(d->getBlock(),type);

		acl::Kernel k(kernelConfig);
		k << i2p.initPosition;
		k << acl::assignmentSafe(d->getEContainer(), f->value(i2p.position));
		k.setup();
		k.compute();
	}
		
	PFRotationField::PFRotationField(acl::VectorOfElements rotationAxis,
	                                 acl::VectorOfElements center):
		axis(rotationAxis),
		c(center)
	{
	}

	acl::VectorOfElements PFRotationField::value(acl::VectorOfElements & pos)
	{
		return crossProduct(axis, pos-c);
	}
		
	SPPositionFunction generatePFRotationField(const AVec<> & a, const AVec<> & c)
	{
		return SPPositionFunction(new PFRotationField(acl::generateVEConstant(a), 
		                                              acl::generateVEConstant(c)));
	}
		
}// asl
