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


#ifndef ACLGENERATORS_H
#define ACLGENERATORS_H

#include "aclMath/aclVectorOfElementsDef.h"
#include "aclMath/aclMatrixOfElements.h"
#include <math/aslVectors.h>
#include <math/aslMatrices.h>


namespace acl
{
	/// Generates VectorOfElements with 1 Element acl::Constant with value a \ingroup generateVE
	template <typename T> VectorOfElements generateVEConstant(T a);
	/// Generates VectorOfElements with 2 Elements acl::Constant with values a and b \ingroup generateVE
	template <typename T> VectorOfElements generateVEConstant(T a, T b);
	/// Generates VectorOfElements with 3 Elements acl::Constant with values a,b and c \ingroup generateVE
	template <typename T> VectorOfElements generateVEConstant(T a, T b, T c);
	/// Generates VectorOfElements with n Elements acl::Constant with values a \ingroup generateVE
	template <typename T> VectorOfElements generateVEConstantN(unsigned int n, T a);
	/// Generates VectorOfElements with \p n Elements acl::Constant with values a[i] \ingroup generateVE
	template <typename T> VectorOfElements generateVEConstant(unsigned int n, const T* const a);
	/// Generates VectorOfElements with a.size() Elements acl::Constant with values a[i] \ingroup generateVE
	template <typename T> VectorOfElements generateVEConstant(const std::vector<T> & a);
	/// Generates VectorOfElements correspondinng to \p a \ingroup generateVE
	template <typename T> VectorOfElements generateVEConstant(const asl::AVec<T> & a);
	/// Generates VectorOfElements correspondinng to \p a \ingroup generateME
	template <typename T> MatrixOfElements generateMEConstant(const asl::AMatr<T> & a);
		
	/// Generates VectorOfElements with \p nComponents Elements acl::Vector with size \p length \ingroup generateVE
	template <typename T> VectorOfElementsData generateVEData(unsigned int length,
	                                                          unsigned int nComponents,
	                                                          CommandQueue queue);

	/// Generates VectorOfElements with \p nComponents Elements acl::Vector with size \p length and default queue \ingroup generateVE
	template <typename T> VectorOfElementsData generateVEData(unsigned int length,
	                                                          unsigned int nComponents = 1);
	
	/// Generates VectorOfElementsData with \p nComponents Elements acl::Array of type \p with size \p length \ingroup generateVE
	VectorOfElementsData generateVEData(unsigned int length,
    	                                TypeID typeID,
		                                unsigned int nComponents,
		                                CommandQueue queue);

		/// Generates VectorOfElementsData with \p nComponents Elements acl::Array of type \p with size \p length and default queue \ingroup generateVE
	VectorOfElementsData generateVEData(unsigned int length,
    	                                TypeID typeID,
		                                unsigned int nComponents = 1);

	
	/// Generates VectorOfElements with \p size Elements acl::LocalArray of type \p typeID with size \p componentSize \ingroup generateVE
	VectorOfElements generateVELocalArray(unsigned int componentSize,
	                                      TypeID typeID,
	                                      unsigned int size);

	/// Generates VectorOfElements with \p size Elements acl::PrivateArray of type \p with data defined by \p data \ingroup generateVE
	template <typename T> VectorOfElements generateVEPrivateArray(const vector<T> & data);
	/// Generates VectorOfElements with \p size Elements acl::PrivateArray of type \p with data defined by \p data \ingroup generateVE
	template <typename T> VectorOfElements generateVEPrivateArray(const vector<asl::AVec<T>> & data);

	/// Generates VectorOfElements with \p size Elements acl::PrivateArray of type \p with size \p componentSize \ingroup generateVE
	template <typename T> VectorOfElements generateVEPrivateArray(const vector<T> & data,
	                                  	                          TypeID typeID);
	/// Generates VectorOfElements with \p size Elements acl::PrivateArray of type \p with size \p componentSize \ingroup generateVE
	template <typename T> VectorOfElements generateVEPrivateArray(const vector<asl::AVec<T>> & data,
	                                  	                          TypeID typeID);
	
	/// Generates VectorOfElements with \p nComponents Elements acl::Subvector with size \p sublength. \p length is the vector size \ingroup generateVE
	template <typename T> VectorOfElements generateVEDataSub(T,
	                                                     unsigned int sublength,
	                                                     unsigned int length,
	                                                     unsigned int nComponents,
	                                                     CommandQueue queue);

	/// Generates VectorOfElements with 1 Element acl::VariableReference with reference on \p a	\ingroup generateVE
	template <typename T> VectorOfElements generateVEVariableR(T& a);
	/// Generates VectorOfElements with 2 Element acl::VariableReference with references on \p a and \p b	\ingroup generateVE
	template <typename T> VectorOfElements generateVEVariableR(T& a, T& b);
	/// Generates VectorOfElements with 3 Element acl::VariableReference with references on \p a, \p b and \p c \ingroup generateVE	
	template <typename T> VectorOfElements generateVEVariableR(T& a, T& b, T& c);
	/// Generates VectorOfElements with nD(a) Element acl::VariableReference with reference on \p a[i] \ingroup generateVE	
	template <typename T> VectorOfElements generateVEVariableR(asl::AVec<T>& a);


	/// Generates VectorOfElements with 1 Element acl::VariableReference with reference on \p a	\ingroup generateVE
	template <typename T> VectorOfElements generateVEVariableSP(std::shared_ptr<T> a);
	/// Generates VectorOfElements with 2 Element acl::VariableReference with references on \p a and \p b	\ingroup generateVE
	template <typename T> VectorOfElements generateVEVariableSP(std::shared_ptr<T> a, std::shared_ptr<T> b);
	/// Generates VectorOfElements with 3 Element acl::VariableReference with references on \p a, \p b and \p c \ingroup generateVE	
	template <typename T> VectorOfElements generateVEVariableSP(std::shared_ptr<T> a, std::shared_ptr<T> b, std::shared_ptr<T> c);
	/// Generates VectorOfElements with nD(a) Element acl::VariableReference with reference on \p a[i] \ingroup generateVE	
	template <typename T> VectorOfElements generateVEVariableSP(std::shared_ptr<asl::AVec<T>> a);
	
	
	/// Generates VectorOfElements with \p n Element of acl::PrivateVariable \ingroup generateVE
	template <typename T> VectorOfElements generateVEPrivateVariable(unsigned int n);

	/// Generates VectorOfElements with \p n Element of acl::PrivateVariable with type \p t \ingroup generateVE
	VectorOfElements generateVEPrivateVariable(unsigned int n, TypeID t);

	/// Generates VectorOfElements with \p n Element of acl::PrivateVariable with type \p t \ingroup generateME
	MatrixOfElements generateMEPrivateVariable(unsigned int nR, unsigned int nC, TypeID t);

	
	/// Generates VectorOfElements which contains SubElements of the corresponding element of \p a \ingroup generateVE
	VectorOfElements generateVESubElements(VectorOfElements a, unsigned int length, int offset);
	
	/// Generates VectorOfElements which contains SubElements of the corresponding element of \p a \ingroup generateVE
	VectorOfElements generateVESubElements(VectorOfElements a, unsigned int length, VectorOfElements offset);

	/// Generates VectorOfElements which contains SubElements of the corresponding element of \p a \ingroup generateVE
	VectorOfElements generateVEShftedElements(VectorOfElements a, int offset);

	/// Generates VectorOfElements which contains SubElements of the corresponding element of \p a \ingroup generateVE
	VectorOfElements generateVEShftedElements(VectorOfElements a, const std::vector<int> & offset);
	
	/// Generates VectorOfElements which contains SubElements of the corresponding element of \p a \ingroup generateVE
	VectorOfElements generateVEShiftedElements(VectorOfElements a, VectorOfElements offset);

	/// Generates VectorOfElements which 1 element correspond to polynom for \p x with \p coefs \ingroup generateVE
	/// The polinom is \f$ x^{n-1}*coefs_0+x^{n-2}coefs_1+...+coefs_{n-1} \f$
	/**
		The polynom contains mad fanction with type specification the type is 
		 defined by \p x  
	*/
	VectorOfElements generateVEPolynom(VectorOfElements x, VectorOfElements coefs);
	
	/// 
	VectorOfElements generateParsedVE(const VectorOfElements & fields,
	                                  const vector<string> & names,
	                                  const string & statement);
		

	/// \ingroup generateME	
	template <typename T=int>MatrixOfElements generateMEUnit(unsigned int n);

	/// \ingroup generateME	
	MatrixOfElements generateMEDiagonal(const VectorOfElements & d);

	/// \ingroup generateME	
	/**
		\param  sc contains sin and cos values
	 */
	MatrixOfElements generateMEGivensRotation(unsigned int k,
	                                          unsigned int l,
	                                          const VectorOfElements & sc);

	
	// insert local variable which takes values accoding to index \ingroup generateVE
	template <typename T> VectorOfElements indexDependedConstant(vector<unsigned int> r,
	                                                             vector<T> values);


	/// Generates VectorOfElements with one Element of type Index
	/// \ingroup generateVE
	VectorOfElements generateVEIndex(unsigned int size = 0);

	/// Generates VectorOfElements with one Element of type GroupID
	/// \ingroup generateVE
	VectorOfElements generateVEGroupID();

	/// Generates VectorOfElements with one Element of type Index
	/// \ingroup generateVE
	VectorOfElements generateVEIndexExt(unsigned int size = 0);
	
}  //namespace acl

#endif // ACLGENERATORS_H
