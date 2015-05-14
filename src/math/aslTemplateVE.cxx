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


#include "aslTemplateVE.h"
#include "acl/acl.h"
#include "acl/aclGenerators.h"
#include <aslUtilities.h>

namespace asl
{
	TemplateVE::TemplateVE()
	{
	}


	void TemplateVE::init(AbstractDataWithGhostNodes & data,
	                      const VectorTemplate & vectorT,
	                      unsigned int k,
	                      bool bIni)
	{
		vectorTemplate=&vectorT;
		vto=vtObject(vectorTemplate);
		unsigned int nV(vectorTemplate->vectors.size());
		
		if (data.getEContainer().size()<1) 
			errorMessage ("TemplateVE: The input data has wrong number of components: " + 
			              numToStr(data.getEContainer().size()));
		if (data.getEContainer().size()<=k) 
			errorMessage ("TemplateVE: The input data has less components than asked: " + 
			              numToStr(data.getEContainer().size())+", k="+ numToStr(k));
		acl::TypeID t(acl::getElementType(data.getEContainer()));

		if (bIni)
		{
			if((values.size()!=nV) || initValues.size()==0)
				copy(acl::generateVEPrivateVariable(nV, t), values);
			initValues.resize(nV);

			using namespace	acl::elementOperators;
			acl::Element e(acl::generateSubElement(data.getEContainer()[k],
			                                       data.getSubContainerSize(),
			                                       data.getSubContainerOffset()));
			for (unsigned int i(0); i < values.size(); ++i)
			{
				int offset(data.getBlock().c2i(vectorTemplate->vectors[i]));
				initValues[i] =generateShiftedElement(e,offset);
			}
			copy(values=initValues,initValues);
		}
		else
		{
			initValues.resize(0);
			values.resize(nV);

			using namespace	acl::elementOperators;
			acl::Element e(acl::generateSubElement(data.getEContainer()[k],
			                                       data.getSubContainerSize(),
			                                       data.getSubContainerOffset()));
			for (unsigned int i(0); i < values.size(); ++i)
			{
				int offset(data.getBlock().c2i(vectorTemplate->vectors[i]));
				values[i] =generateShiftedElement(e,offset);
			}
		}
	}

	void TemplateVE::init(DistanceFunction & data,
	                      acl::VectorOfElements & position,
	                      const VectorTemplate & vectorT,
	                      unsigned int k)
	{
		vectorTemplate=&vectorT;
		vto=vtObject(vectorTemplate);
		unsigned int nV(vectorTemplate->vectors.size());
		initValues.resize(nV);
		acl::TypeID t(acl::getElementType(position));

		if(values.size()!=nV)
			copy(acl::generateVEPrivateVariable(nV, t), values);

		using namespace	acl::elementOperators;
		for (unsigned int i(0); i < values.size(); ++i)
			initValues[i] =data.getDistance(position + vectorTemplate->vectors[i])[0];
		copy(values=initValues,initValues);
	}
	
	TemplateVE::TemplateVE(AbstractDataWithGhostNodes & data, 
	                       const VectorTemplate & vectorT, 
	                       unsigned int k,
	                       bool bIni)
	{
		init(data,vectorT,k,bIni);
	}

	TemplateVE::TemplateVE(DistanceFunction & data,
	                       acl::VectorOfElements & position, 
	                       const VectorTemplate & vectorT, 
	                       unsigned int k)
	{
		init(data, position,vectorT,k);
	}
	
	TemplateVE::TemplateVE(const acl::VectorOfElements & val, 
	                       const VectorTemplate & vectorT):
		vectorTemplate(&vectorT),
		vto(vtObject(&vectorT)),
		values(val),
		initValues(acl::generateVEConstantN (vectorTemplate->vectors.size(),0.))
	{
	}

		
	acl::VectorOfElements TemplateVE::getValue(unsigned int i)
	{
		return subVE(values,i);
	}

	TemplateVE operator+ (const TemplateVE &a, const TemplateVE &b)
	{
		if (a.vectorTemplate != b.vectorTemplate)
			errorMessage ("Two TemplateVE corespond to different templates");
		return TemplateVE (a.values+b.values,*a.vectorTemplate);
	}
	TemplateVE operator- (const TemplateVE &a, const TemplateVE &b)	
	{
		if (a.vectorTemplate != b.vectorTemplate)
			errorMessage ("Two TemplateVE corespond to different templates");
		return TemplateVE (a.values-b.values,*a.vectorTemplate);
	}
	TemplateVE operator* (const TemplateVE &a, const TemplateVE &b)	
	{
		if (a.vectorTemplate != b.vectorTemplate)
			errorMessage ("Two TemplateVE corespond to different templates");		
		return TemplateVE (productOfElements (a.values,b.values),*a.vectorTemplate);
	}

	TemplateVE operator/ (const TemplateVE &a, const TemplateVE &b)
	{
		if (a.vectorTemplate != b.vectorTemplate)
			errorMessage ("Two TemplateVE corespond to different templates");		
		return TemplateVE (divisionOfElements (a.values,b.values),*a.vectorTemplate);
	}
		           
	acl::VectorOfElements laplas(const TemplateVE & a)
	{
		using namespace acl;
		return a.values*generateVEConstant (a.vectorTemplate->laplasCoefs);
	}

	acl::VectorOfElements dx(const TemplateVE & a)
	{
		using namespace acl;
		unsigned int nD(a.vectorTemplate->numberOfDimentions());
		if (nD<1)
			errorMessage ("(dx) The dimensionality lees than 1");		

		return a.values*generateVEConstant(a.vectorTemplate->dxCoefs[0]);
	}

	acl::VectorOfElements dy(const TemplateVE & a)
	{
		using namespace acl;
		unsigned int nD(a.vectorTemplate->numberOfDimentions());
		if (nD<2)
			errorMessage ("(dy) The dimensionality lees than 2");		

		return a.values*generateVEConstant(a.vectorTemplate->dxCoefs[1]);
	}

	acl::VectorOfElements dz(const TemplateVE & a)
	{
		using namespace acl;
		unsigned int nD(a.vectorTemplate->numberOfDimentions());
		if (nD<3)
			errorMessage ("(dz) The dimensionality lees than 3");		

		return a.values*generateVEConstant(a.vectorTemplate->dxCoefs[2]);
	}
		
	acl::VectorOfElements dxi(const TemplateVE & a, unsigned int i)
	{
		using namespace acl;
		unsigned int nD(a.vectorTemplate->numberOfDimentions());
		if (nD<i)
			errorMessage ("(dxi) The dimensionality lees than "+ numToStr(i));		

		return a.values*generateVEConstant(a.vectorTemplate->dxCoefs[i]);
	}
		
	acl::VectorOfElements gradient(const TemplateVE & a)
	{
		using namespace acl;
		VectorOfElements r(a.values * 
		                   generateVEConstant (a.vectorTemplate->dxCoefs[0]));
		unsigned int nD(a.vectorTemplate->numberOfDimentions());
		for (unsigned int i(1); i < nD; ++i)
			copy(cat(r,
			         a.values * 
			         generateVEConstant(a.vectorTemplate->dxCoefs[i])),
			     r);
		return r;
	}

	acl::VectorOfElements div(const TemplateVE & ax, 
	                          const TemplateVE & ay)
	{
		return dx(ax)+dy(ay);
	}

	acl::VectorOfElements div(const TemplateVE & ax, 
	                          const TemplateVE & ay, 
	                          const TemplateVE & az)
	{
		return dx(ax)+dy(ay)+dz(az);
	}

	acl::VectorOfElements div(const vector<TemplateVE> & a)
	{
		using namespace acl;
		unsigned int nD(a[0].vectorTemplate->numberOfDimentions());
		if (nD!=a.size())
			errorMessage ("(div) The dimensionality of the underline template does not equal to the number of components");		

		auto res(dxi(a[0], 0));
		for(unsigned int i(1); i<nD; ++i)         
			copy(res + dxi(a[i], i), res);
		return res;
	}
		
	acl::VectorOfElements divProduct(const vector<TemplateVE> & a, const TemplateVE & c)
	{
		unsigned int nD(a[0].vectorTemplate->numberOfDimentions());
		unsigned int nv(a[0].values.size());
		
//		return a.values*generateVEConstant(a.vectorTemplate->dxCoefs[i]);
		vector<acl::VectorOfElements> a0(nD);
		for(unsigned int i(0); i<nD; ++i)
			copy(catN(subVE(a[i].values,0),nv), a0[i]);
		acl::VectorOfElements c0(catN(subVE(c.values,0),nv));
		auto res(0.5 * productOfElements((a[0].values + a0[0]), (c.values + c0)) * 
		         acl::generateVEConstant(a[0].vectorTemplate->dxCoefs[0]));
		for(unsigned int i(1); i<nD; ++i)         
			copy(res + 
			     0.5 * productOfElements((a[i].values + a0[i]), (c.values + c0)) * 
		         acl::generateVEConstant(a[0].vectorTemplate->dxCoefs[i]), res);
		return res;
	}
		
	acl::VectorOfElements divAgradB(const TemplateVE & a, const TemplateVE & b)
	{
		if (a.vectorTemplate != b.vectorTemplate) 
			errorMessage ("divAgradB: \"a\" and \"b\" have different templates");

		unsigned int nv(a.values.size());
		acl::VectorOfElements a0(catN(subVE(a.values,0),nv));
		acl::VectorOfElements b0(catN(subVE(b.values,0),nv));
		acl::VectorOfElements coefs(acl::generateVEConstant (a.vectorTemplate->gradientCoefs));
		
		return productOfElements(coefs,a.values+a0)*(b.values-b0);
	}

	acl::VectorOfElements dIdJ(unsigned int i, unsigned int j, const TemplateVE & a)
	{
		using namespace acl;
		acl::VectorOfElements a0(catN(subVE(a.values,0),a.values.size()));
		acl::VectorOfElements coefs(acl::generateVEConstant (a.vectorTemplate->dIdJCoefs[i][j]));
		acl::VectorOfElements coefLap(acl::generateVEConstant (a.vectorTemplate->dIdJLapCoef));
		if (i==j)
			return coefs*(a.values-a0)-coefLap*laplas(a);
		return coefs*(a.values-a0);
	}
		
	acl::VectorOfElements interpolate(const TemplateVE & a, acl::VectorOfElements e)
	{
		auto & vt(*a.vectorTemplate); 
		unsigned int nd(nD(vt));
		AVec<> one(nd,1.);
		acl::VectorOfElements res(1);
		copy(productOfElements((one - AVec<>(vt.vectors[0])) + 
		                        acl::productOfElements(acl::generateVEConstant(2.* AVec<>(vt.vectors[0]) -one), e)) * 
		     subVE(a.values,0),
		     res); 
		for(unsigned int i(1); i < a.values.size(); ++i)
			copy(productOfElements((one - AVec<>(vt.vectors[i])) + 
			                        acl::productOfElements(acl::generateVEConstant(2.* AVec<>(vt.vectors[i]) -one), e)) * 
				 subVE(a.values,i) + 
			     res,
			 	 res); 
			
		return res;
	}

		
}// asl
