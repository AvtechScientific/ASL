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


#include "aclGenerators.h"
#include "aclMath/aclVectorOfElements.h"
#include "acl.h"
#include "DataTypes/aclConstant.h"
#include "DataTypes/aclArray.h"
#include "DataTypes/aclSubvector.h"
#include "DataTypes/aclVariableReference.h"
#include "DataTypes/aclVariableSP.h"
#include "DataTypes/aclPrivateVariable.h"
#include "DataTypes/aclPrivateArray.h"
#include "DataTypes/aclIndex.h"
#include "DataTypes/aclIndexExt.h"
#include "DataTypes/aclGroupID.h"

#include "aclTypesList.h"


namespace acl
{

	template <typename T> VectorOfElements generateVEConstant(T a)
	{
		VectorOfElements v(1); 
		v[0] = Element(new Constant<T>(a));
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEConstant(t a);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEConstant(T a, T b)
	{
		VectorOfElements v(2); 
		v[0] = Element(new Constant<T>(a));
		v[1] = Element(new Constant<T>(b));
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEConstant(t a, t b);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEConstant(T a, T b, T c)
	{
		VectorOfElements v(3); 
		v[0] = Element(new Constant<T>(a));
		v[1] = Element(new Constant<T>(b));
		v[2] = Element(new Constant<T>(c));
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEConstant(t a, t b, t c);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression


	template <typename T> VectorOfElements generateVEConstantN(unsigned int n, T a)
	{
		VectorOfElements v(n);
		v[0] = Element(new Constant<T>(a));
		for (unsigned int i(1); i < n; ++i)
			v[i] = v[0];
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEConstantN(unsigned int n, t a);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEConstant(unsigned int n,
	                                                          const T* const a)
	{
		VectorOfElements v(n);
		for (unsigned int i(0); i < n; ++i)
			v[i] = Element(new Constant<T>(a[i]));
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEConstant(unsigned int n, const t* const a);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEConstant(const vector<T> & a)
	{
		return generateVEConstant(a.size(), &a[0]);		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEConstant(const vector<t> & a);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEConstant(const asl::AVec<T> &a)
	{
		VectorOfElements v(a.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i)
			v[i] = Element(new Constant<T>(a[i]));
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEConstant(const asl::AVec<t> &a);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> MatrixOfElements generateMEConstant(const asl::AMatr<T> &a)
	{
		MatrixOfElements m(a.getNRow(), a.getNCol());
		copy(generateVEConstant(a.getInternalVec()),
		     m.getInternalVector());
		return m;		
	}

	template MatrixOfElements generateMEConstant(const asl::AMatr<double> &a);
	template MatrixOfElements generateMEConstant(const asl::AMatr<float> &a);
		
	template <typename T> 
	VectorOfElementsData generateVEData(unsigned int length,
	                                    unsigned int nComponents,
	                                    CommandQueue queue)
	{
		VectorOfElementsData v(nComponents);
		for (unsigned int i(0); i < nComponents; ++i)
			v[i] = ElementData(new Array<T>(length, queue));
		return v;		
	}	

	#define BOOST_TT_rep_expression(r, data, t) \
	template VectorOfElementsData generateVEData<t>(unsigned int length, \
	                                                unsigned int nComponents, \
	                                                CommandQueue queue);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> 
	VectorOfElementsData generateVEData(unsigned int length,
	                                    unsigned int nComponents)
	{
		return generateVEData<T>(length,nComponents,hardware.defaultQueue);
	}

	#define BOOST_TT_rep_expression(r, data, t) \
	template VectorOfElementsData generateVEData<t>(unsigned int length, \
	                                                unsigned int nComponents);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
	
	VectorOfElementsData generateVEData(unsigned int length,
    	                                TypeID typeID,
		                                unsigned int nComponents,
		                                CommandQueue queue)
	{
		VectorOfElementsData v(nComponents);
		for (unsigned int i(0); i < nComponents; ++i)
			v[i] = generateElementArray(typeID, length, queue);
		return v;
	}
	

	VectorOfElements generateVELocalArray(unsigned int componentSize,
	                                      TypeID typeID,
	                                      unsigned int size)
	{
		VectorOfElements v(size);
		for (unsigned int i(0); i < size; ++i)
			v[i] = generateElementLocalArray(typeID, componentSize);
		return v;
	}

	template <typename T> VectorOfElements generateVEPrivateArray(const vector<T> & data)
	{
		VectorOfElements v(1);
		v[0] = Element(new PrivateArray<T>(data));
		return v;		
	}
	#define BOOST_TT_rep_expression(r, data, T) \
	template VectorOfElements generateVEPrivateArray(const vector<T> & d);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEPrivateArray(const vector<asl::AVec<T>> & data)
	{
		unsigned int nd(data[0].getSize());
		VectorOfElements v(nd);
		vector<T> d(data.size());
		for(unsigned int i(0); i < nd; ++i)
		{
			for(unsigned int j(0); j < data.size(); ++j)
				d[j] = data[j][i];
			v[i] = Element(new PrivateArray<T>(d));
		}
		return v;		
	}
	#define BOOST_TT_rep_expression(r, data, T) \
	template VectorOfElements generateVEPrivateArray(const vector<asl::AVec<T>> & d);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

		
	template <typename T> VectorOfElements generateVEPrivateArray(const vector<T> & d,
	                                  	                          TypeID typeID)
	{
		VectorOfElements v(1);
		switch (typeID) 
		{
			#define BOOST_TT_rep_expression(r, data, tt) \
				case typeToTypeID<tt>():	\
					{   \
						vector<tt> newD(d.size());  \
						copy(d.begin(),d.end(),newD.begin());   \
						copy(generateVEPrivateArray(newD), v);  \
					}   \
					break;	
			BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
			#undef BOOST_TT_rep_expression	
		}
		return v;
	}
	#define BOOST_TT_rep_expression(r, data, T) \
	template VectorOfElements generateVEPrivateArray(const vector<T> & d, TypeID typeID);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
	
	template <typename T> VectorOfElements generateVEPrivateArray(const vector<asl::AVec<T>> & d,
	                                  	                          TypeID typeID)
	{
		VectorOfElements v(1);
		switch (typeID) 
		{
			#define BOOST_TT_rep_expression(r, data, tt) \
				case typeToTypeID<tt>():	\
					{   \
						vector<asl::AVec<tt>> newD(d.size());  \
						copy(d.begin(),d.end(),newD.begin());   \
						copy(generateVEPrivateArray(newD), v);  \
					}   \
					break;	
			BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
			#undef BOOST_TT_rep_expression	
		}
		return v;
	}
	#define BOOST_TT_rep_expression(r, data, T) \
	template VectorOfElements generateVEPrivateArray(const vector<asl::AVec<T>> & d, TypeID typeID);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression


		
	template <typename T> 
	VectorOfElements generateVEDataSub(T,
	                                   unsigned int sublength,
	                                   unsigned int length,
	                                   unsigned int nComponents,
	                                   CommandQueue queue)
	{
		VectorOfElements v(nComponents);
		for (unsigned int i(0); i < nComponents; ++i)
		{
			shared_ptr<Array<T> > vec(new Array<T>(length, queue));
			v[i] = Element(new Subvector<T>(vec, sublength, 0u));
		}
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
	template VectorOfElements generateVEDataSub(t,	\
	                                            unsigned int sublength,	\
	                                            unsigned int length,	\
	                                            unsigned int nComponents,	\
	                                            CommandQueue queue);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEVariableR(T& a)
	{
		VectorOfElements v(1); 
		v[0] = Element(new VariableReference<T>(a));
		return v;
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEVariableR(t & a);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEVariableR(T& a, T& b)
	{
		VectorOfElements v(2); 
		v[0] = Element(new VariableReference<T>(a));
		v[1] = Element(new VariableReference<T>(b));
		return v;	
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEVariableR(t & a, t & b);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
	
	template <typename T> VectorOfElements generateVEVariableR(T& a, T& b, T& c)
	{
		VectorOfElements v(3); 
		v[0] = Element(new VariableReference<T>(a));
		v[1] = Element(new VariableReference<T>(b));
		v[2] = Element(new VariableReference<T>(c));
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEVariableR(t & a, t & b, t & c);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEVariableR(asl::AVec<T>& a)
	{
		unsigned int n(nD(a));
		VectorOfElements v(n); 
		for (unsigned int i(0); i < n; ++i)
			v[i] = Element(new VariableReference<T>(a[i]));
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEVariableR(asl::AVec<t> & a);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEVariableSP(std::shared_ptr<T> a)
	{
		VectorOfElements v(1); 
		v[0] = Element(new VariableSP<T>(a));
		return v;
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEVariableSP(std::shared_ptr<t> a);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> VectorOfElements generateVEVariableSP(std::shared_ptr<T> a,
	                                                            std::shared_ptr<T> b)
	{
		VectorOfElements v(2); 
		v[0] = Element(new VariableSP<T>(a));
		v[1] = Element(new VariableSP<T>(b));
		return v;
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEVariableSP(std::shared_ptr<t> a, std::shared_ptr<t> b);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression	

	template <typename T> VectorOfElements generateVEVariableSP(std::shared_ptr<T> a,std::shared_ptr<T> b, std::shared_ptr<T> c)
	{
		VectorOfElements v(3); 
		v[0] = Element(new VariableSP<T>(a));
		v[1] = Element(new VariableSP<T>(b));
		v[2] = Element(new VariableSP<T>(c));
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEVariableSP(std::shared_ptr<t> a, std::shared_ptr<t> b, std::shared_ptr<t> c);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression	

	template <typename T> VectorOfElements generateVEVariableSP(std::shared_ptr<asl::AVec<T>> a)
	{
		unsigned int n(nD(*a));
		VectorOfElements v(n); 
		for (unsigned int i(0); i < n; ++i)
			v[i] = Element(new VariableSP<T>(std::shared_ptr<T>(a, & (*a)[i]) ));
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEVariableSP(std::shared_ptr<asl::AVec<t>> a);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression	

	
	template <typename T> VectorOfElements generateVEPrivateVariable(unsigned int n)
	{
		VectorOfElements v(n); 
		for (unsigned int i(0); i < n; ++i)
			v[i] = Element(new PrivateVariable<T>());
		return v;		
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template VectorOfElements generateVEPrivateVariable<t>(unsigned int i);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression	

	VectorOfElements generateVEPrivateVariable(unsigned int n, TypeID t)
	{
		VectorOfElements v;
		switch (t) 
		{
			#define BOOST_TT_rep_expression(r, data, tt) \
				case typeToTypeID<tt>():	\
					copy(generateVEPrivateVariable<tt>(n), v);	\
					break;	
			BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
			#undef BOOST_TT_rep_expression	
		}
		return v;		
	}

	MatrixOfElements generateMEPrivateVariable(unsigned int nR, unsigned int nC, TypeID t)
	{
		MatrixOfElements res(nR,nC);
		copy(generateVEPrivateVariable(nR*nC, t), res.getInternalVector());
		return res;
	}

	
	VectorOfElements generateVESubElements(VectorOfElements a,
	                                       unsigned int length,
	                                       int offset)
	{
		unsigned int n(a.size());
		VectorOfElements v(n); 
		for (unsigned int i(0); i < n; ++i)
			v[i] = generateSubElement(a[i], length, offset);
		return v;		
	}


	VectorOfElements generateVESubElements(VectorOfElements a,
	                                       unsigned int length,
	                                       VectorOfElements offset)
	{
		unsigned int n(a.size());
		if (n != offset.size())
			errorMessage("generateVESubElements(): two VectorOfElements have different sizes");
		VectorOfElements v(n);
		for (unsigned int i(0); i < n; ++i)
			v[i] = generateSubElement(a[i], length, offset[i]);
		return v;		
	}


	VectorOfElements generateVEShiftedElements(VectorOfElements a,
	                                          int offset)
	{
		unsigned int n(a.size());
		VectorOfElements v(n); 
		for (unsigned int i(0); i < n; ++i)
			v[i] = generateShiftedElement(a[i], offset);
		return v;		
	}

	VectorOfElements generateVEShiftedElements(VectorOfElements a,
	                                           VectorOfElements offset)
	{
		unsigned int n(a.size());
		if (n != offset.size())
			errorMessage("generateVEShiftedElements(): two VectorOfElements have different sizes");
		VectorOfElements v(n); 
		for (unsigned int i(0); i < n; ++i)
			v[i] = generateShiftedElement(a[i], offset[i]);
		return v;		
	}


	VectorOfElements generateVEShiftedElements(VectorOfElements a,
	                                           const vector<int> & offset)
	{
		unsigned int n(a.size());
		if (n != offset.size())
			errorMessage("generateVEShiftedElements(): VectorOfElements and vector<int> have different sizes");
		VectorOfElements v(n); 
		for (unsigned int i(0); i < n; ++i)
			v[i] = generateShiftedElement(a[i], offset[i]);
		return v;		
	}


	VectorOfElements generateParsedVE(const VectorOfElements & fields,
	                                  const vector<string> & names,
	                                  const string & statement)
	{
		if (fields.size() != names.size())
			errorMessage("generateParsedVE(): VectorOfElements \"fields\" and vector<string> \"names\" have different sizes");
		vector<pair<Element, string> > f(names.size());
		for (unsigned int i(0); i < names.size(); ++i)
		{
			f[i].first = fields[i];
			f[i].second = names[i];
		}

		VectorOfElements a(1u);
		a[0] = elementOperators::parse(f, statement);
		return a;			
	}
		
	VectorOfElements generateVEPolynom(VectorOfElements x,
	                                   VectorOfElements coefs)
	{
		if (x.size() > 1)
			errorMessage("generateVEPolynom(): \"x\" has size more than 1");
		if (coefs.size() < 1)
			errorMessage("generateVEPolynom(): size of \"coefs\" less than 1");
		auto type(getElementType(x));
		VectorOfElements p(x.size());
		copy(subVE(coefs, 0), p);
		for (unsigned int i(1); i < coefs.size(); ++i)
			copy(mad(p, x, subVE(coefs, i), type), p);
//			p[0] = p[0] * x[0] + coefs[i];
		return p;	
	}

	template <typename T> MatrixOfElements generateMEUnit(unsigned int n)
	{
		MatrixOfElements m(n,n);
		Element c0(new Constant<T>(0));
		Element c1(new Constant<T>(1));
		
		for(unsigned int i(0); i<n; ++i)
			for(unsigned int j(0); j<n; ++j)
				m.setElement(i,j,c0);
		for(unsigned int i(0); i<n; ++i)
			m.setElement(i,i,c1);
		return m;
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template MatrixOfElements generateMEUnit<t>(unsigned int n);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression	
	
	MatrixOfElements generateMEDiagonal(const VectorOfElements & d)
	{
		unsigned int n(d.size());
		MatrixOfElements m(n,n);

		Element c0(new Constant<int>(0));
		
		for(unsigned int i(0); i<n; ++i)
			for(unsigned int j(0); j<n; ++j)
				m.setElement(i,j,c0);
		for(unsigned int i(0); i<n; ++i)
			m.setElement(i,i,d[i]);
		
		return m;
	}
	
/*	template <typename T> VectorOfElements indexDependedConstant(vector<unsigned int> r, vector<T> values)
	{
		Index in;
		
		Element e(select())
		return 
	}*/


	VectorOfElements generateVEIndex(unsigned int size)
	{
		VectorOfElements v(1);
		Element ind(new Index(size));
		v[0] = ind;
		return v;
	}

	VectorOfElements generateVEGroupID()
	{
		VectorOfElements v(1);
		v[0].reset(new GroupID());
		return v;
	}

	VectorOfElements generateVEIndexExt(unsigned int size)
	{
		VectorOfElements v(1);
		Element ind(new IndexExt(size));
		v[0] = ind;
		return v;
	}
	
} // acl
