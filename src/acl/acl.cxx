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


#include "acl.h"
#include "Operators/aclElementAssignmentSafe.h"
#include "Operators/aclElementIfElse.h"
#include "Operators/aclElementFor.h"
#include "Operators/aclElementSum.h"
#include "Operators/aclElementSubtraction.h"
#include "Operators/aclElementProduct.h"
#include "Operators/aclElementDivision.h"
#include "Operators/aclElementSin.h"
#include "Operators/aclElementCos.h"
#include "Operators/aclElementSqrt.h"
#include "Operators/aclOperatorGeneric.h"
#include "Operators/aclElementGenericBinary.h"
#include "Operators/aclGenericAtomicFunction.h"
#include "Operators/aclElementGenericUnary.h"
#include "Operators/aclElementSelect.h"
#include "Operators/aclElementMad.h"
#include "Operators/aclElementExcerpt.h"
#include "Operators/aclElementParser.h"
#include "Operators/aclElementSyncCopy.h"
#include "Operators/aclElementConvert.h"
#include "../aslUtilities.h"
#include "Kernels/aclKernel.h"
#include "DataTypes/aclConstant.h"
#include "DataTypes/aclIndex.h"
#include "DataTypes/aclVariableReference.h"
#include "DataTypes/aclArray.h"
#include "DataTypes/aclLocalArray.h"
#include "Kernels/aclExpressionContainer.h"
#include "aclMath/aclVectorOfElementsDef.h"
#include "aclMath/aclMatrixOfElements.h"

#include<string>

using namespace std;
using namespace asl;

namespace acl
{
	namespace elementOperators
	{

		Element operator-(Element e)
		{
			return Element(new ElementGenericUnary(e, "-", false));
		}
		
		Element operator+(Element e1, Element e2)
		{
			return Element(new ElementSum(e1, e2));
		}


		Element operator-(Element e1, Element e2)
		{
			return Element(new ElementSubtraction(e1, e2));
		}


		Element operator*(Element e1, Element e2)
		{
			return Element(new ElementProduct(e1, e2));
		}


		Element operator/(Element e1, Element e2)
		{
			return Element(new ElementDivision(e1, e2));
		}

		
		Element operator%(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "%"));
		}

		
		Element operatorAssignmentSafe(Element e1, Element e2)
		{
			return Element(new ElementAssignmentSafe(e1, e2));
		}


		Element operatorAssignment(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "="));
		}

		
		Element operator+=(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "+="));
		}


		Element operator-=(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "-="));
		}


		Element operator*=(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "*="));
		}


		Element operator/=(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "/="));
		}


		Element operator>(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, ">"));
		}


		Element operator<(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "<"));
		}

	
		Element operator>=(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, ">="));
		}

	
		Element operator<=(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "<="));
		}

	
		Element isEqual(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "=="));
		}

	
		Element isNotEqual(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, " != "));
		}


		Element operator||(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "||"));
		}


		Element operator&&(Element e1, Element e2)
		{
			return Element(new ElementGenericBinary(e1, e2, "&&"));
		}

		Element operator!(Element e)
		{
			return Element(new ElementGenericUnary(e, "!", false));
		}

		Element convert(const TypeID tName, Element e1, bool strong)
		{
			return Element(new ElementConvert(e1, tName, strong));
		}
		

		Element select(Element e1, Element e2, Element e3)
		{
			return Element(new ElementSelect(e1,e2,e3));
		}
	
		Element mad(Element e1, Element e2, Element e3)
		{
			return Element(new ElementMad(e1, e2, e3));
		}

		Element sin(Element e)
		{
			return Element(new ElementSin(e));
		}


		Element cos(Element e)
		{
			return Element(new ElementCos(e));
		}


		Element sqrt(Element e)
		{
			return Element(new ElementSqrt(e));
		}

		Element rsqrt(Element e)
		{
			return Element(new ElementGenericUnary(e, "native_rsqrt"));
		}
	

		Element log(Element e)
		{
			return Element(new ElementGenericUnary(e, "native_log"));
		}


		Element log10(Element e)
		{
			return Element(new ElementGenericUnary(e, "log10"));
		}


		Element powI(Element a, unsigned int i)
		{
			Element n(new Constant<int>(i));
			return Element(new ElementGenericBinaryFunction(a, n, "pown"));
		}


		Element exp(Element a)
		{
			return Element(new ElementGenericUnary(a, "exp"));
		}


		Element fabs(Element a)
		{
			return Element(new ElementGenericUnary(a, "fabs"));
		}

		Element abs(Element a)
		{
			return Element(new ElementGenericUnary(a, "abs"));
		}

		Element abs_diff(Element a, Element b)
		{
			return Element(new ElementGenericBinaryFunction(a, b, "abs_diff"));
		}
		
		Element floor(Element a)
		{
			return Element(new ElementGenericUnary(a, "floor"));
		}
		
		Element isnan(Element a)
		{
			return Element(new ElementGenericUnary(a, "isnan"));
		}

		Element copysign(Element a, Element b)
		{
			return Element(new ElementGenericBinaryFunction(a, b, "copysign"));
		}

		Element sign(Element a)
		{
			return Element(new ElementGenericUnary(a, "sign"));
		}
		
		Element min(Element a, Element b)
		{
			return Element(new ElementGenericBinaryFunction(a, b, "min"));
		}

		Element max(Element a, Element b)
		{
			return Element(new ElementGenericBinaryFunction(a, b, "max"));
		}
		
		Element atomic_add(Element e1, Element e2)
		{
			return Element(new ElementGenericAtomicFunction(e1, e2, "atomic_add"));
		}


		Element atomic_sub(Element e1, Element e2)
		{
			return Element(new ElementGenericAtomicFunction(e1, e2, "atomic_sub"));
		}


		Element atomic_xchg(Element e1, Element e2)
		{
			return Element(new ElementGenericAtomicFunction(e1, e2, "atomic_xchg"));
		}

			
		Element excerpt(Element source, Element filter)
		{
			return Element(new ElementExcerpt(source, filter));
		}
		

		Element parse(const vector<pair<Element, string> > & elementNamePairs, const string & statement)
		{
			shared_ptr<ElementParser> parser(new ElementParser);
			parser->setStatement(statement);

			for (unsigned int i = 0; i < elementNamePairs.size(); ++i)
			{
				parser->addElementNamePair(elementNamePairs[i].first, elementNamePairs[i].second);
			}
			return parser;
		}


		Element printfFunction(string args)
		{
			return Element(new OperatorGeneric("printf(" + args + ")"));
		}

		
		Element syncCopy(Element source,
		                 Element destination,
		                 Element srcOffset,
		                 Element dstOffset,
		                 Element length)
		{
			return Element(new ElementSyncCopy(source, destination,
			                                   srcOffset, dstOffset,
			                                   length));
		}

		Element any(Element e)
		{
			return Element(new ElementGenericUnarySIMD(e, "any"));
		}

		Element all(Element e)
		{
			return Element(new ElementGenericUnarySIMD(e, "all"));
		}
		

		/// Sets work-group barrier
		Element barrier(string flags)
		{
			return Element(new OperatorGeneric("barrier(" + flags + ")"));
		}

		
		Element ifElse(Element condition, 
		               const vector<Element> & thenBody, 
		               const vector<Element> & elseBody)
		{
			std::shared_ptr<ElementIfElse> elIf(new ElementIfElse(condition));
			for(unsigned int i(0); i < thenBody.size(); ++i)
				elIf->addBodyExpressionIf(thenBody[i]);
			for(unsigned int i(0); i < elseBody.size(); ++i)
				elIf->addBodyExpressionIf(elseBody[i]);

			return elIf;
		}


		Element forLoop(Element initialization,
		                Element condition,
		                Element increase,
		                const vector<Element> & body)
		{
			std::shared_ptr<ElementFor> elFor(new ElementFor (initialization, condition, increase));
			for(unsigned int i(0); i < body.size(); ++i)
			{
				elFor->addBodyExpression(body[i]);
			}
			return Element(elFor);
		}


		Element returnStatement()
		{
			return Element(new OperatorGeneric("return"));
		}

		Element nan(TypeID t)
		{
			std::string s; 
			switch (t)
			{
				case TYPE_DOUBLE:
					s=std::string("nan(ulong(1))");
					break;	
				case TYPE_FLOAT:
					s="nan(uint(1))";
					break;
				default:
					errorMessage("nan: the input variable has an uncorrect type");
			}
			return Element(new OperatorGeneric(s));
		}

	}

	
	//	RTTI functions

	bool isConstant()
	{
		return false;
	}
	

	bool isSingleValue()
	{
		return false;
	}
	

	bool isMemBlock(Element e)
	{
		return (dynamic_cast<MemBlock*>(e.get()) != NULL);
	}

	/// Copies source to destination, resizes destination to accommodate source.
	template <typename T> void copy(MemBlock & source, T* destination)
	{
		cl_int status = 0;
		cl::Event event;
		// Enqueue readBuffer. Blocking read (CL_TRUE)!
		status = source.getQueue()->enqueueReadBuffer(source.getBuffer(),
		                                              CL_TRUE,
		                                              0,
		                                              source.getSize() * sizeof(T),
		                                              destination,
		                                              NULL,
		                                              &event);
		
		errorMessage(status, "queue::enqueueReadBuffer()");
		status = event.wait();
		errorMessage(status, "Event::wait() - event");
	}
	
	/// Copies source to destination, resizes destination to accommodate source.
	template <typename T> void copy(MemBlock & source, std::vector<T> & destination)
	{
		// this line is necesary in order to avoid unnecesary memory copy trigered by resize()
		destination.clear();
		destination.resize(source.getSize());
		copy(source, &(destination[0]));
	}

	template <typename T> void copy(Element source, std::vector<T> & destination)
	{
		if (isMemBlock(source))
			copy(dynamic_cast<MemBlock& >(*source), destination);
		else
			errorMessage("copy() failed. First argument is not a MemBlock or has unproper type");
	}

	template void copy(Element source, std::vector<cl_int> &destination);
	template void copy(Element source, std::vector<cl_float> &destination);
	template void copy(Element source, std::vector<cl_double> &destination);

	template <typename T> void copy(Element source, T* destination)
	{
		if (isMemBlock(source))
			copy(dynamic_cast<MemBlock& >(*source),destination);
		else
			errorMessage("copy() failed. First argument is not a MemBlock or has unproper type");
	}

	template void copy(Element source, cl_int* destination);
	template void copy(Element source, cl_float* destination);
	template void copy(Element source, cl_double* destination);
	
	
	template <typename T> void copy(T* source, MemBlock & destination)
	{
		cl_int status = 0;
		cl::Event event;		
		// Write data to buffer 
		status = destination.getQueue()->enqueueWriteBuffer(destination.getBuffer(),
		                                                    CL_TRUE,
		                                                    0,
		                                                    destination.getSize() * sizeof(T),
		                                                    &source[0],
		                                                    NULL,
		                                                    &event);
		errorMessage(status, "copy() - queue::enqueueWriteBuffer()");
		status = event.wait();
		errorMessage(status, "Event::wait() - event");
		
	}

	
	template <typename T> void copy(std::vector<T> & source, MemBlock & destination)
	{
		if (source.size() == destination.getSize())
			copy(&source[0], destination);
		else
			errorMessage("copy() - write to MemBlock failed. Sizes do not match");
	}


	template <typename T> void copy(MemBlock & source, MemBlock & destination)
	{
		if (source.size() == destination.getSize())
		{
			if (source.getQueue().get() == destination.getQueue().get())
			{
				// \todoIgal swap buffers???
			}
			else
			{
				shared_ptr<void> buffer = source.map();
				copy((T *)buffer, destination);
			}
		}
		else
		{
			errorMessage("copy() - write from MemBlock to MemBlock failed. Sizes do not match");
		}
		    
	}
	
	template <typename T> void copy(std::vector<T> & source, Element destination)
	{
		if (isMemBlock(destination))
			copy(source, dynamic_cast<MemBlock& >(*destination));				
		else
			errorMessage("copy() - Second argument is not a MemBlock type or has unproper type");
	}

	template void copy(std::vector<cl_int> & source, Element destination);
	template void copy(std::vector<cl_float> & source, Element destination);
	template void copy(std::vector<cl_double> & source, Element destination);

	template <typename T> void copy(T* source, Element destination)
	{
		if (isMemBlock(destination))
			copy(source, dynamic_cast<MemBlock& >(*destination));				
		else
			errorMessage("copy() - Second argument is not a MemBlock type or has unproper type");
	}

	template void copy(cl_int* source, Element destination);
	template void copy(cl_uint* source, Element destination);
	template void copy(cl_float* source, Element destination);
	template void copy(cl_double* source, Element destination);
	template void copy(cl_long* source, Element destination);

	std::vector<Element> & operator<<(std::vector<Element> & ec,
	                                 const std::vector<Element> & a)
	{
		for (unsigned int i(0); i < a.size(); ++i)
			ec.push_back(a[i]);
		return ec;
	}

	ExpressionContainer & operator<<(ExpressionContainer & ec,
	                                 const std::vector<Element> & a)
	{
		for (unsigned int i(0); i < a.size(); ++i)
			ec.addExpression(a[i]);
		return ec;
	}
	
	ExpressionContainer & operator<<(ExpressionContainer & ec,
	                                 const MatrixOfElements & a)
	{
		ec << (vector<Element>) (a.getInternalVector()) ;
		return ec;
	}
	
	ExpressionContainer & operator<<(ExpressionContainer & ec,
	                                 const ExpressionContainer & a)
	{
		ec << a.expression;
		return ec;
	}
	
	void initData(Element a,
	              Element initializationValue,
	              const KernelConfiguration & kernelConfig)
	{
		Kernel k(kernelConfig);
		{ 
			using namespace elementOperators;
			k.addExpression(operatorAssignment(a, initializationValue));
		}
		k.compute();
	}


	Element generateSubElement(Element e, unsigned int size, int offset)
	{
		Element ind(new Index(size));
		Element offs(new Constant<cl_int>(offset));
		using namespace elementOperators;
		Element res(new ElementExcerpt(e, ind + offs));
		return res;
	}


	Element generateSubElement(Element e, unsigned int size, int * offset)
	{
		Element ind(new Index(size));
		Element offs(new VariableReference<cl_int>(*offset));
		using namespace elementOperators;
		Element res(new ElementExcerpt(e, ind + offs));
		return res;
	}


	Element generateSubElement(Element e, unsigned int size, Element offset)
	{
		Element ind(new Index(size));
		using namespace elementOperators;
		Element res(new ElementExcerpt(e, ind + offset));
		return res;
	}	


	Element generateShiftedElement(Element e,int offset)
	{
		Element ind(new Index(e->getSize()));
		Element offs(new Constant<cl_int>(offset));
		using namespace elementOperators;
		Element res(new ElementExcerpt(e, ind + offs));
		return res;
	}


	Element generateShiftedElement(Element e, unsigned int size, int * offset)
	{
		Element ind(new Index(e->getSize()));
		Element offs(new VariableReference<cl_int>(*offset));
		using namespace elementOperators;
		Element res(new ElementExcerpt(e, ind + offs));
		return res;
	}


	Element generateShiftedElement(Element e, Element offset)
	{
		Element ind(new Index(e->getSize()));
		using namespace elementOperators;
		Element res(new ElementExcerpt(e, ind + offset));
		return res;
	}


	ElementData generateElementArray(TypeID typeID,
	                                 unsigned int size,
	                                 CommandQueue q)
	{
		ElementData v;
		switch (typeID) 
		{
			case TYPE_DOUBLE:
				v.reset(new Array<cl_double>(size, q));
				break;	
			case TYPE_FLOAT:
				v.reset(new Array<cl_float>(size, q));
				break;
		    case TYPE_INT:
				v.reset(new Array<cl_int>(size, q));		
				break;
		    case TYPE_UINT:
				v.reset(new Array<cl_uint>(size, q));		
				break;
		    case TYPE_LONG:
				v.reset(new Array<cl_long>(size, q));		
		}
		return v;
	}

	ElementData generateElementArray(TypeID typeID,
	                                 unsigned int size)
	{
		return generateElementArray(typeID, size, hardware.defaultQueue);
	}

	Element generateElementLocalArray(TypeID typeID,
	                                  unsigned int size)
	{
		Element v;
		switch (typeID) 
		{
			case TYPE_DOUBLE:
				v.reset(new LocalArray<cl_double>(size));
				break;	
			case TYPE_FLOAT:
				v.reset(new LocalArray<cl_float>(size));
				break;
		    case TYPE_INT:
				v.reset(new LocalArray<cl_int>(size));		
				break;
		    case TYPE_UINT:
				v.reset(new LocalArray<cl_uint>(size));		
				break;
		    case TYPE_LONG:
				v.reset(new LocalArray<cl_long>(size));		
		}
		return v;
	}	

} // namespace acl
