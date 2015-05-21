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


#ifndef ACLVECTORELEMENTSOPERATIONS_H
#define ACLVECTORELEMENTSOPERATIONS_H

#include "aclVectorOfElementsDef.h"

namespace acl
{
	
	class ElementBase;
	class MemBlock;
	class KernelConfiguration;


	typedef std::shared_ptr<ElementBase> Element;
	typedef std::shared_ptr<MemBlock> ElementData;

	class VectorOfElements;

	typedef shared_ptr<VectorOfElementsData> SPVectorOfElementsData;
	typedef shared_ptr<VectorOfElements> SPVectorOfElements;

	/// Creates VectorOfElementsData with same structure ElementData objects as \p a
	/// \relates VectorOfElementsData
	VectorOfElementsData clone(VectorOfElementsData a);

	/// Creates VectorOfElementsData with same structure ElementData objects as \p a and with \p n numbers of elements
	/// \relates VectorOfElementsData
	VectorOfElementsData clone(VectorOfElementsData a, unsigned int n);
	
	
	/// Writes data from \p initializationValue to \p a
	/// \relates VectorOfElements
	void initData(VectorOfElements a,
	              VectorOfElements initializationValue,
	              const KernelConfiguration & kernelConfig);

	/// Writes data from \p initializationValue to \p a
	/// \relates VectorOfElements
	void initData(VectorOfElements a,
	              VectorOfElements initializationValue);
	

	/// Copies the VectorOfElements class. 
	/** \relates VectorOfElements 
		the \p destination is resized automaticaly
	*/
	void copy(const vector<Element> & source,
	          VectorOfElements & destination);

	/// Copies the VectorOfElements class. 
	/** \relates VectorOfElements 
		the \p destination is resized automaticaly
	*/
	void copy(const vector<ElementData> & source,
	          VectorOfElements & destination);
	
	/// Copies the VectorOfElements class. 
	/** \relates VectorOfElements 
		the \p destination is resized automaticaly
		 \p start and \p end define elements to copy: [start: end]
	*/
	void copy(const vector<Element> & source,
	          VectorOfElements & destination, unsigned int start, unsigned int end);
	

	/// Copies the VectorOfElementsData class. 
	/** \relates VectorOfElementsData
		the \p destination is resized automaticaly
	*/
	void copy(const VectorOfElementsData & source,
	          VectorOfElementsData & destination);


	/// returns VectorOfElements class containing several elements of \p a. 
	/// \relates VectorOfElements 
	VectorOfElements subVE(const VectorOfElements & source,
	                       unsigned int start,
	                       unsigned int end);

	/// returns VectorOfElements class containing one element \p i of \p a. 
	/// \relates VectorOfElements 
	VectorOfElements subVE(const VectorOfElements & source,
	                       unsigned int i);

	/// returns VectorOfElements class containing elements with numbers \p iList of \p source. 
	/// \relates VectorOfElements 
	VectorOfElements subVE(const VectorOfElements & source,
	                       const vector<unsigned int> & iList);
	
	/// returns VectorOfElements class containing several elements of \p a. 
	/// \relates VectorOfElements 
	VectorOfElementsData subVE(const VectorOfElementsData & source,
	                           unsigned int start,
	                           unsigned int end);


	/// Swaps buffers between corresponding
	/// elements of two VectorOfElementsData classes. 
	/// \relates VectorOfElementsData
	void swapBuffers(const VectorOfElementsData & a,
	                 const VectorOfElementsData & b);

	///\defgroup VectorOfElementsOp VectorOfElements Operations 

	///\defgroup VectorOfElementsVectorOfElementsOp  VectorOfElements-VectorOfElements Operations 	
	///@{ 
		
	/// Creates VectorOfElements containing assignment operation result \p a=\p b  
	/// \relates VectorOfElements
	VectorOfElements assignmentSafe(const VectorOfElements & a,
	                                const VectorOfElements & b);	
	

	/// Creates VectorOfElements containing operation result of element \p a  
	///	\relates VectorOfElements
	VectorOfElements operator-(const VectorOfElements & a);	
	

	/// Creates VectorOfElements containing
	/// operation result of elements of \p a and \p b 
	///	\relates VectorOfElements
	VectorOfElements operator+=(const VectorOfElements & a,
	                            const VectorOfElements & b);


	/// Creates VectorOfElements containing operation result of elements of \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator-=(const VectorOfElements & a,
	                            const VectorOfElements & b);


	/// Creates VectorOfElements containing operation result of elements of \p a and \p b 
	/// \relates VectorOfElements
	///	It is assumed that size of \p b is 1
	VectorOfElements operator*=(const VectorOfElements & a,
	                            const VectorOfElements & b);


	/// Creates  VectorOfElements containing operation result of elements of \p a and \p b 
	/// \relates VectorOfElements
	/// It is supposed that size of \p b is 1
	VectorOfElements operator/=(const VectorOfElements & a,
	                            const VectorOfElements & b);


	/// Creates VectorOfElements containing operation result of elements of \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator+(const VectorOfElements & a,
	                           const VectorOfElements & b);


	/// Creates VectorOfElements containing operation result of elements of \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator-(const VectorOfElements & a,
	                           const VectorOfElements & b);


	/// Creates  VectorOfElements containing operation result of elements of \p a and \p b 
	/// \relates VectorOfElements
	/// In case of sizes of \p a and \p b are equal the product is a scalar (dot)
	/// product. If one of sizes is 1 than the function results product of a scalar on vector
	VectorOfElements operator*(const VectorOfElements & a,
	                           const VectorOfElements & b);


	/// Creates VectorOfElements containing operation result of elements of \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator/(const VectorOfElements & a,
	                           const VectorOfElements & b);

	/// Creates VectorOfElements containing operation result of elements of \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator%(const VectorOfElements & a,
	                           const VectorOfElements & b);
	

	/// Creates VectorOfElements containing operation result of elements \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator==(const VectorOfElements & a,
	                            const VectorOfElements & b);


	/// Creates VectorOfElements containing operation result of elements \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator!=(const VectorOfElements & a,
	                            const VectorOfElements & b);


	/// Creates VectorOfElements containing operation result of elements  \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator>(const VectorOfElements & a, const VectorOfElements & b);


	/// Creates VectorOfElements containing operation result of elements \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator<(const VectorOfElements & a,
	                           const VectorOfElements & b);

	/// Creates VectorOfElements containing operation result of elements \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator<=(const VectorOfElements & a,
	                            const VectorOfElements & b);

	/// Creates VectorOfElements containing operation result of elements \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements operator>=(const VectorOfElements & a,
	                            const VectorOfElements & b);
	
	/// Creates VectorOfElements containing operation result of elements \p a && \p b 
	/// \relates VectorOfElements
	VectorOfElements operator&&(const VectorOfElements & a,
	                           const VectorOfElements & b);

	/// Creates VectorOfElements containing operation result of elements \p a || \p b 
	/// \relates VectorOfElements
	VectorOfElements operator||(const VectorOfElements & a,
	                           const VectorOfElements & b);

	/// Creates VectorOfElements containing operation result of element \p a  
	///	\relates VectorOfElements
	VectorOfElements operator!(const VectorOfElements & a);	

	
	/// Creates VectorOfElements \p c which elements are crossproduct of corresponding elements of \p a and \p b 
	/// \relates VectorOfElements
	/// The function is defined for 2D and 3D cases
	VectorOfElements crossProduct(const VectorOfElements & a,
	                              const VectorOfElements & b);


	/// Creates VectorOfElements \p b corresponding to a scala product a*a 
	/// \relates VectorOfElements
	/// The function is defined for 2D and 3D cases
	inline VectorOfElements l2(const VectorOfElements & a);
	

	/// Creates VectorOfElements \p c which
	/// elements are product of corresponding elements of \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements productOfElements(const VectorOfElements & a,
	                                   const VectorOfElements & b);

	/// Creates VectorOfElements \p c which
	/// element is a product of all elements of \p a 
	/// \relates VectorOfElements
	VectorOfElements productOfElements(const VectorOfElements & a);

	/// Creates VectorOfElements \p c which elements are division of corresponding elements of \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements divisionOfElements(const VectorOfElements & a,
	                                    const VectorOfElements & b);
	
	/// Creates VectorOfElements \p c which elements are min of corresponding elements of \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements min(const VectorOfElements & a,
	                     const VectorOfElements & b);

	/// Creates VectorOfElements \p c which elements are min of corresponding elements of \p a and \p b with type specification 
	/// \relates VectorOfElements
	VectorOfElements min(const VectorOfElements & a,
	                     const VectorOfElements & b,
	                     TypeID type);
	
	/// Creates VectorOfElements \p c which elements are min of corresponding elements of \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements minAbs(const VectorOfElements & a,
	                        const VectorOfElements & b);
	
	/// Creates VectorOfElements \p c which elements are max of corresponding elements of \p a and \p b 
	/// \relates VectorOfElements
	VectorOfElements max(const VectorOfElements & a,
	                     const VectorOfElements & b);

	/// Creates VectorOfElements \p c which elements are max of corresponding elements of \p a and \p b with type specification 
	/// \relates VectorOfElements
	VectorOfElements max(const VectorOfElements & a,
	                     const VectorOfElements & b,
	                     TypeID type);

	///returns a with sign of b copysign(a,b)
	/// \relates VectorOfElements
	VectorOfElements copysign(const VectorOfElements & a,
	                          const VectorOfElements & b);

	///returns a with sign of b copysign(a,b) with type specification
	/// \relates VectorOfElements
	VectorOfElements copysign(const VectorOfElements & a,
	                          const VectorOfElements & b,
	                          TypeID t);

	///returns sign of a 
	/// \relates VectorOfElements
	VectorOfElements sign(const VectorOfElements & a);
	
	/// Creates VectorOfElements which elements are result of excerpt
	/// operation of \p source and \p filter elements
	/// \relates VectorOfElements
	/// The \p filter shoule have only 1 component
	VectorOfElements excerpt(const VectorOfElements & source,
	                         const VectorOfElements & filter);
		 
	
	/// Creates VectorOfElements containing operation
	/// result of elements \p a, \p b and \p c 
	/// For a scalar type, result = c ? b : a
	/** 
		 \relates VectorOfElements
		 The function can operate in two modes: \p c has 1 components and 
		 several components. 
	*/
	VectorOfElements select(const VectorOfElements & a,
	                        const VectorOfElements & b,
	                        const VectorOfElements & c);

	/// Creates VectorOfElements containing operation
	/// result of elements \p a, \p b and \p c 
	/// For a scalar type, result = c ? b : a, with type specification
	/** 
		 \relates VectorOfElements
		 The function can operate in two modes: \p c has 1 components and 
		 several components. 
	*/
	VectorOfElements select(const VectorOfElements & a,
	                        const VectorOfElements & b,
	                        const VectorOfElements & c, 
	                        TypeID t);

	/// Creates VectorOfElements containing operation
	/// result of elements \p a, \p b and \p c 
	/// For a scalar type, result = c ? b : 0, with type specification
	/** 
		 \relates VectorOfElements
		 The function can operate in two modes: \p c has 1 components and 
		 several components. 
	*/
	VectorOfElements select(const VectorOfElements & b,
	                        const VectorOfElements & c, 
	                        TypeID t);
	
	/// Creates VectorOfElements containing operation
	/// result of elements \p a, \p b and \p c 
	/// For a scalar type, result = a * b + c
	/// \relates VectorOfElements
	VectorOfElements mad(const VectorOfElements & a,
	                     const VectorOfElements & b,
	                     const VectorOfElements & c);
	
	/// Creates VectorOfElements containing operation
	/// result of elements \p a, \p b and \p c with type specification 
	/// For a scalar type, result = a * b + c
	/// \relates VectorOfElements
	VectorOfElements mad(const VectorOfElements & a,
	                     const VectorOfElements & b,
	                     const VectorOfElements & c,
	                     TypeID t);
	
	/// Creates VectorOfElements containing result log(a); a should contain only one element 
	/// \relates VectorOfElements
	VectorOfElements log(const VectorOfElements & a);


	/// Creates VectorOfElements containing result log10(a); a should contain only one element 
	/// \relates VectorOfElements
	VectorOfElements log10(const VectorOfElements & a);

	
	/// Creates VectorOfElements containing result a^i; a should contain only one element 
	/// \relates VectorOfElements
	VectorOfElements powI(const VectorOfElements & a, unsigned int i);


	/// Creates VectorOfElements containing result exp(a); a should contain only one element 
	/// \relates VectorOfElements
	VectorOfElements exp(const VectorOfElements & a);

	/// Creates VectorOfElements containing result exp(a); a should contain only one element 
	/// \relates VectorOfElements
	VectorOfElements sqrt(const VectorOfElements & a);

	/// Creates VectorOfElements containing result exp(a); a should contain only one element 
	/// \relates VectorOfElements
	VectorOfElements rsqrt(const VectorOfElements & a);
	
	/// Creates VectorOfElements containing result element-wise operation fabs(a);
	/// \relates VectorOfElements
	VectorOfElements fabs(const VectorOfElements & a);

	/// Creates VectorOfElements containing result element-wise operation abs(a);
	/// \relates VectorOfElements
	VectorOfElements abs(const VectorOfElements & a);

	/// Creates VectorOfElements containing result element-wise operation abs_diff(a, b);
	/// \relates VectorOfElements
	VectorOfElements abs_diff(const VectorOfElements & a, const VectorOfElements & b);
	
	/// Creates VectorOfElements containing result element-wise operation \f$ floor(a_i) \f$ 
	/// \relates VectorOfElements
	VectorOfElements floor(const VectorOfElements & a);
	
	/// Creates VectorOfElements with openCL type casting;
	/// \relates VectorOfElements
	VectorOfElements convert(acl::TypeID type, const VectorOfElements & a, bool strong=true);
	
	/// Creates VectorOfElements containing a min element 
	/// \relates VectorOfElements
	VectorOfElements minElement(const VectorOfElements & a);

	/// Creates VectorOfElements containing an element with minimal absolute value 
	/// \relates VectorOfElements
	VectorOfElements minAbsElement(const VectorOfElements & a);
	
	/// Creates VectorOfElements containing a min element 
	/// \relates VectorOfElements
	VectorOfElements maxElement(const VectorOfElements & a);

	/// Creates VectorOfElements containing a sum of elements 
	/// \relates VectorOfElements
	VectorOfElements sumOfElements(const VectorOfElements & a);

	/// Creates VectorOfElements containing a && operation result of elements 
	/// \relates VectorOfElements
	VectorOfElements andOfElements(const VectorOfElements & a);

	/// Creates VectorOfElements containing a || operation result of elements 
	/// \relates VectorOfElements
	VectorOfElements orOfElements(const VectorOfElements & a);
	
	/// Concatinates two vectors \p a and \p b
	/// \relates VectorOfElements
	/**	\f$ cat\left(
		          \left[\begin{array}{c}
		             a_1\\a_2\\...\\a_n\\
		          \end{array}\right],
		          \left[\begin{array}{c}
		             b_1\\b_2\\...\\b_m\\
		          \end{array}\right]
		          \right)=
		 		          \left[\begin{array}{c}
		             a_1\\a_2\\ \vdots \\a_n\\
	                 b_1\\b_2\\ \vdots \\b_m\\
		          \end{array}\right]. \f$		 
	 */
	VectorOfElements cat(const VectorOfElements & a,
	                     const VectorOfElements & b);

	/// Concatinates two vectors \p a and \p b
	/// \relates VectorOfElements
	/**	\f$ cat\left(
		          \left[\begin{array}{c}
		             a_1\\a_2\\...\\a_n\\
		          \end{array}\right],
		          \left[\begin{array}{c}
		             b_1\\b_2\\...\\b_m\\
		          \end{array}\right]
		          \right)=
		 		          \left[\begin{array}{c}
		             a_1\\a_2\\ \vdots \\a_n\\
	                 b_1\\b_2\\ \vdots \\b_m\\
		          \end{array}\right]. \f$		 
	 */
	VectorOfElementsData cat(const VectorOfElementsData & a,
	                         const VectorOfElementsData & b);

	/// Concatinates three vectors \p a  \p b and \p c
	/// \relates VectorOfElements
	/**	\f$ cat\left(
		          \left[\begin{array}{c}
		             a_1\\...\\a_n\\
		          \end{array}\right],
		          \left[\begin{array}{c}
		             b_1\\...\\b_m\\
		          \end{array}\right],
		          \left[\begin{array}{c}
		             c_1\\...\\c_l\\
		          \end{array}\right]
		          \right)=
		 		          \left[\begin{array}{c}
		             a_1\\...\\a_n\\
	                 b_1\\...\\b_m\\
	                 c_1\\...\\c_l\\
		          \end{array}\right]. \f$		 
		 
	 */
	VectorOfElements cat(const VectorOfElements & a,
	                     const VectorOfElements & b,
	                     const VectorOfElements & c);

	/// Concatinates \p n VectorOfElements \p a   
	/// \relates VectorOfElements
	VectorOfElements cat(const VectorOfElements * a, 
	                     unsigned int n);

	/// Concatinates VectorOfElements \p a to itself \p n times
	/// \relates VectorOfElements
	VectorOfElements catN(const VectorOfElements & a, 
	                      unsigned int n);

	
	/// Creates VectorOfElements containing operation result of elements of \p a and a double 
	///	\relates VectorOfElements
	/// The function creates automaticaly a constant of a corresponding type
	template <typename T> VectorOfElements operator+=(const VectorOfElements & a,
	                                                  const T & b);


	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	///	\relates VectorOfElements
	/// The function creates automaticaly a constant of a corresponding type
	template <typename T> VectorOfElements operator-=(const VectorOfElements & a,
	                                                  const T & b);


	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator*=(const VectorOfElements & a,
	                                                  const T & b);


	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator/=(const VectorOfElements & a,
	                                                  const T & b);
	

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator+(const VectorOfElements & a,
	                                                 const T & b);
	

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator+(const T & a,
	                                                 const VectorOfElements & b);
	

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	///	\relates VectorOfElements
	/// The function creates automaticaly a constant of a corresponding type
	template <typename T> VectorOfElements operator-(const VectorOfElements & a,
	                                                 const T & b);
	

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	///	\relates VectorOfElements
	/// The function creates automaticaly a constant of a corresponding type
	template <typename T> VectorOfElements operator-(const T & a,
	                                                 const VectorOfElements & b);


	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	///	\relates VectorOfElements
	/// The function creates automaticaly a constant of a corresponding type
	template <typename T> VectorOfElements operator*(const VectorOfElements & a,
	                                                 const T & b);
	

	/// Creates VectorOfElements containing
	/// operation result of elements of \p a and a double 
	///	\relates VectorOfElements
	/// The function creates automaticaly a constant of a corresponding type
	template <typename T> VectorOfElements operator*(const T & a,
	                                                 const VectorOfElements & b);


	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator/(const VectorOfElements & a,
	                                                 const T & b);
	
	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator%(const VectorOfElements & a,
	                                                 const T & b);

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator/(const T & b,
	                                                 const VectorOfElements & a);
	
	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator%(const T & b,
	                                                 const VectorOfElements & a);

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator>(const VectorOfElements & a,
	                                                 const T & b);
	

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator>(const T & b,
	                                                 const VectorOfElements & a);
	

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator<(const VectorOfElements & a,
	                                                 const T & b);
	

	/// Creates VectorOfElements containing
	/// operation result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator<(const T & b,
	                                                 const VectorOfElements & a);

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator>=(const VectorOfElements & a,
	                                                  const T & b);
	

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator>=(const T & b,
	                                                  const VectorOfElements & a);

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator<=(const VectorOfElements & a,
	                                                  const T & b);
	

	/// Creates VectorOfElements containing operation
	/// result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator<=(const T & b,
	                                                  const VectorOfElements & a);

	/// Creates VectorOfElements containing
	/// operation result of elements of \p a and a double 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator==(const VectorOfElements & a,
	                                                  const T & b);
	

	/// Creates VectorOfElements containing
	/// operation result of elements of \p a and a double \p b 
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator==(const T & b,
	                                                  const VectorOfElements & a);
	

	/// Creates VectorOfElements containing
	/// operation result of elements of \p a and a double \p b
	/**	\relates VectorOfElements
		 The function creates automaticaly a constant of a corresponding type
	*/
	template <typename T> VectorOfElements operator!=(const VectorOfElements & a,
	                                                  const T & b);
	

	/// Creates VectorOfElements containing
	/// operation result of elements of \p a and a double \p b  
	///	\relates VectorOfElements
	/// The function creates automaticaly a constant of a corresponding type
	template <typename T> VectorOfElements operator!=(const T & b,
	                                                  const VectorOfElements & a);		
	///@} 

	/// generates code corresponding \f$ \vec a /= |\vec a|\f$ 	
	vector<Element> gcNormalize(const VectorOfElements & a);
	/// generates code corresponding \f$ |\vec a|^2 \f$. Result will be stored in \p l2 	
	vector<Element> gcLength2(const VectorOfElements & a, const VectorOfElements & l2);
	/// generates code corresponding \f$ |\vec a| \f$. Result will be stored in \p l 	
	vector<Element> gcLength(const VectorOfElements & a, const VectorOfElements & l);
	
		
//------------------------------ Implementation ----------------

	inline VectorOfElements l2(const VectorOfElements & a)
	{
		return a * a;
	}

		
}  //namespace acl

#endif // ACLVECTORELEMENTSOPERATIONS_H
