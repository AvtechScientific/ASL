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


#ifndef ACL_H
#define ACL_H

#include "aclStdIncludes.h"
#include <memory>
#include "Kernels/aclKernelConfiguration.h"
//#include "aclHardware.h"
#include "aclTypes.h"

namespace cl
{
	class CommandQueue;
}


/// Advanced Computational Language
namespace acl
{
	class ExpressionContainer;
	class VectorOfElements;
	class MatrixOfElements;
	class MemBlock;
	class ElementBase;
	template <class T> class Array;
	typedef std::shared_ptr<MemBlock> ElementData;
	typedef std::shared_ptr<ElementBase> Element;
	extern const KernelConfiguration KERNEL_BASIC;
	typedef std::shared_ptr<cl::CommandQueue> CommandQueue;	

	
	/// definitions of mathematical operators and functions for class ::Element
	namespace elementOperators
	{
	// Math Operators
		/**
			 \ingroup MathOperators
		 */
		Element operator-(Element e);

		/**
			 \ingroup MathOperators
		 */
		Element operator+(Element e1, Element e2);

		/**
			 \ingroup MathOperators
		 */
		Element operator-(Element e1, Element e2);

		/**
			 \ingroup MathOperators
		 */
		Element operator*(Element e1, Element e2);

		/**
			 \ingroup MathOperators
		 */
		Element operator/(Element e1, Element e2);

		/**
			 \ingroup MathOperators
		 */
		Element operator%(Element e1, Element e2);


		
	// Assignment Operators
		/**
			 \ingroup AssignmentOperators
		 */		
		Element operatorAssignment(Element e1, Element e2);

		/**
			 \ingroup AssignmentOperators
		 */		
		Element operatorAssignmentSafe(Element e1, Element e2);

		/**
			 \ingroup AssignmentOperators
		 */		
		Element operator+=(Element e1, Element e2);

		/**
			 \ingroup AssignmentOperators
		 */		
		Element operator-=(Element e1, Element e2);

		/**
			 \ingroup AssignmentOperators
		 */		
		Element operator*=(Element e1, Element e2);

		/**
			 \ingroup AssignmentOperators
		 */		
		Element operator/=(Element e1, Element e2);


		/**
			 \ingroup BooleanOperators
		 */		
		Element operator>(Element e1, Element e2);

		/**
			 \ingroup BooleanOperators
		 */	
		Element operator<(Element e1, Element e2);

		/**
			 \ingroup BooleanOperators
		 */	
		Element operator>=(Element e1, Element e2);

		/**
			 \ingroup BooleanOperators
		 */	
		Element operator<=(Element e1, Element e2);

		/**
			 \ingroup BooleanOperators
		 */	
		Element isEqual(Element e1, Element e2);

		/**
			 \ingroup BooleanOperators
		 */	
		Element isNotEqual(Element e1, Element e2);

		/**
			 \ingroup BooleanOperators
		 */	
		Element operator&&(Element e1, Element e2);

		/**
			 \ingroup BooleanOperators
		 */	
		Element operator||(Element e1, Element e2);

		/**
			 \ingroup MathOperators
		 */
		Element operator!(Element e);

		
	// Mathematical Functions
		/**
			 \ingroup MathFunctions
		 */
		Element sin(Element e);

		/**
			 \ingroup MathFunctions
		 */
		Element cos(Element e);

		/**
			 \ingroup MathFunctions
		 */
		Element sqrt(Element e);

		/**
			 \ingroup MathFunctions
		 */
		Element rsqrt(Element e);

		/**
			 \ingroup MathFunctions
		 */
		Element log(Element e);

		/**
			 \ingroup MathFunctions
		 */
		Element log10(Element e);
		
		/**
			 a^i
			 \ingroup MathFunctions
		 */
		Element powI(Element a, unsigned int i);

		/**
			 exp(a)
			 \ingroup MathFunctions
		 */
		Element exp(Element a);

		/**
			 fabs(a)
			 \ingroup MathFunctions
		 */
		Element fabs(Element a);

		/**
			 abs(a)
			 \ingroup MathFunctions
		 */
		Element abs(Element a);

		/**
			 abs_diff(a,b)
			 \ingroup MathFunctions
		 */
		Element abs_diff(Element a,Element b);
		
		/**
			 floor(a)
			 \ingroup MathFunctions
		 */
		Element floor(Element a);
		
		/**
			 isnan(a)
			 \ingroup MathFunctions
		 */
		Element isnan(Element a);

		/// Return nan of corresponding float type
		/// \ingroup MathFunctions
		Element nan(TypeID t);
		
		/**
		 	returns \p a with sign of \p b
			 copysign(a,b)
			 \ingroup MathFunctions
		 */
		Element copysign(Element a,Element b);

		/**
		 	returns 1 with sign of \p a and 0 for a = +-0 
			 copysign(a,b)
			 \ingroup MathFunctions
		 */
		Element sign(Element a);
		
		/**
			 min(a,b)
			 \ingroup MathFunctions
		 */
		Element min(Element a,Element b);

		/**
			 max(a,b)
			 \ingroup MathFunctions
		 */
		Element max(Element a,Element b);
		
		/// Atomic sum
		/// \ingroup AtomicFunctions
		Element atomic_add(Element e1, Element e2);

		/// Atomic subtraction
		/// \ingroup AtomicFunctions
		Element atomic_sub(Element e1, Element e2);

		/// Atomic exchange
		/// \ingroup AtomicFunctions
		Element atomic_xchg(Element e1, Element e2);


		/// ternary mad operator e1 * e2 + e3
		/// \ingroup MathFunctions
		Element mad(Element e1, Element e2, Element e3);


	// Synchronization Functions
		/**
		 Synchronous copy: Global memory <-> Local memory

		 \ingroup SynchronizationFunctions
		 */
		Element syncCopy(Element source,
		                 Element destination,
		                 Element srcOffset,
		                 Element dstOffset,
		                 Element length);



		/**
		 Sets work-group barrier

		 \ingroup SynchronizationFunctions
		 */
		Element barrier(std::string flags = "CLK_LOCAL_MEM_FENCE");



	// Special Purpose Functions
		/**
		 Creates excerpt from source defined by filter
		 by replacing of "INDEX" occurrences
		 in source->str() by filter->str()
			 
		 \ingroup SpecialPurposeFunctions
		 */
		Element excerpt(Element source, Element filter);

		/**
		 Parses statement by replacing occurrences of
		 element's names by corresponding element->str()
		 from elementNamePairs
			 
		 \ingroup SpecialPurposeFunctions
		 */
		Element parse(const std::vector<std::pair<Element, std::string> > & elementNamePairs, 
		              const std::string & statement);

		/**
		 type conversion
		 \ingroup SpecialPurposeFunctions
		 */			
		Element convert(const TypeID tName, Element e1, bool strong=true);

		/// printf()
		/// \ingroup SpecialPurposeFunctions
		Element printfFunction(std::string args);
		
		
		/// Return statement (terminates the execution of a kernel)
		/// \ingroup ControlStructures
		Element returnStatement();

		
		/// If-Else conditional structure
		/// \ingroup ControlStructures
		Element ifElse(Element condition, 
		               const std::vector<Element> & thenBody, 
		               const std::vector<Element> & elseBody);

			
		/// ternary branching operator
		/// \ingroup ControlStructures
		Element select(Element e1, Element e2, Element e3);


		/// For loop
		/// \ingroup ControlStructures
		Element forLoop(Element initialization, 
		                Element condition, 
		                Element increase, 
		                const std::vector<Element> & body);

		/// Corresponds to the openCL operation any	
		Element any(Element e);
		/// Corresponds to the openCL operation all
		Element all(Element e);
		

	}

	
	//	RTTI functions

	bool isConstant(Element e);

	/// The function returns true when the input is a single valued object e.g. aclConstatnt, aclVariable 
	bool isSingleValue(Element e); 

	bool isMemBlock(Element e);

	/// Copies source to destination, resizes destination to accommodate source.
	/// \ingroup HostInteractionFunctions
	template <typename T> void copy(MemBlock &source, T* destination);

	/// Copies source to destination, exit(1) if sizes do not match.
	/// \ingroup HostInteractionFunctions
	template <typename T> void copy(T* source, MemBlock &destination);

	/// Copies source to destination, resizes destination to accommodate source.
	/// \ingroup HostInteractionFunctions
	template <typename T> void copy(MemBlock &source, std::vector<T> &destination);

	/// Copies source to destination, exit(1) if sizes do not match.
	/// \ingroup HostInteractionFunctions
	template <typename T> void copy(std::vector<T> &source, MemBlock &destination);

	/// Copies source to destination, exit(1) if sizes do not match.
	/// \ingroup HostInteractionFunctions
	template <typename T> void copy(MemBlock &source, MemBlock &destination);

	/// Copies source to destination, resizes destination to accommodate source.
	/// \ingroup HostInteractionFunctions
	template <typename T> void copy(Element source, std::vector<T> &destination);	
	
	/// Copies source to destination, exit(1) if sizes do not match.
	/// \ingroup HostInteractionFunctions
	template <typename T> void copy(std::vector<T> &source, Element destination);

	/// Copies source to destination.
	/// \ingroup HostInteractionFunctions
	template <typename T> void copy(Element source, T* destination);	
	
	/// Copies source to destination.
	/// \ingroup HostInteractionFunctions
	template <typename T> void copy(T* source, Element destination);

	/// puts a vector<Element> in ExpressionContainer
	/**
		 \relates ExpressionContainer
		 \relates VectorOfElements
	 */
	ExpressionContainer & operator<<(ExpressionContainer & ec,
	                                 const std::vector<Element> & a);
	
	/// puts a vector<Element> in ExpressionContainer
	/**
		 \relates ExpressionContainer
		 \relates VectorOfElements
	 */
	std::vector<Element> & operator<<(std::vector<Element> & ec,
	                                 const std::vector<Element> & a);

	
	/// puts a MatrixOfElements in ExpressionContainer
	/**
		 \relates ExpressionContainer
		 \relates VectorOfElements
	 */
	ExpressionContainer & operator<<(ExpressionContainer & ec,
	                                 const MatrixOfElements & a);
	
	/// puts a ExpressionContainer in ExpressionContainer
	/**
		 \relates ExpressionContainer
		 \relates VectorOfElements
	 */
	ExpressionContainer & operator<<(ExpressionContainer & ec,
	                                 const ExpressionContainer & a);
	
	/// function writes data from \p initializationValue to \p a
	///	\ingroup SpecialPurposeFunctions
	void initData(Element a,
	              Element initializationValue,
	              const KernelConfiguration & kernelConfig = KERNEL_BASIC);

	/// function creates subElement with given length and offset; offset is constant
	Element generateSubElement(Element, unsigned int size, int offset);		
	
	///function creates subElement with given length and offset; offset can be a variable
	Element generateSubElement(Element, unsigned int size, int * offset);		

	Element generateSubElement(Element, unsigned int size, Element offset);		

	Element generateShiftedElement(Element, int offset);			
	Element generateShiftedElement(Element, int* offset);		
	Element generateShiftedElement(Element, Element offset);		

	template <typename T> void swapBuffers(std::shared_ptr<Array<T> >a, 
	                                       std::shared_ptr<Array<T> > b);


	ElementData generateElementArray(TypeID typeID,
	                                 unsigned int size);
	ElementData generateElementArray(TypeID typeID,
	                                 unsigned int size,
	                                 CommandQueue queue_);


	Element generateElementLocalArray(TypeID typeID,
	                                  unsigned int size);

} // namespace acl
#endif // ACL_H
