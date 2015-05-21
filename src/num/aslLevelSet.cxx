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


#include "aslLevelSet.h"
#include "acl/acl.h"
#include "acl/aclHardware.h"
#include "acl/Kernels/aclKernel.h"
#include "acl/Kernels/aclKernelConfigurationTemplates.h"
#include "acl/aclGenerators.h"
#include "acl/aclMath/aclMathAlg.h"
#include <algorithm>
#include <math/aslTemplatesExtras.h>
#include <math/aslTemplateVE.h>

namespace asl
{

	LevelSet::~LevelSet()
	{}
		
	LevelSet::LevelSet():
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		distanceField(),
		distanceFieldInternalData(),
		vectorTemplate(NULL),
		vto(NULL)
	{
	}
	
	LevelSet::LevelSet(Data df):
		kernel(new acl::Kernel(acl::KERNEL_SIMDUA)),
		distanceField(df),
		distanceFieldInternalData(),
//		vectorTemplate(nearestNeigboursVT(nD(df->getBlock())))
		vectorTemplate(nearestNeigboursPVT(nD(df->getBlock()))),
		vto(vtObject(vectorTemplate))
	{		
	}	

	void LevelSet::initKernelPropagation()
	{
		distanceTVE.reset(new TemplateVE(*distanceField, *vectorTemplate));
		(*kernel) << distanceTVE->initValues;

		initVelocityComputation();
		initDistancesComputation();
			
		kernel->setup();
		cout<<"!!!!!!! "<<acl::getKernelPrivateMemSize(*kernel)<<endl;
		cout<<"!!!!!!! "<<acl::getKernelLocalMemSize(*kernel)<<endl;
//		cout<<kernel->getKernelSource()<<endl;
		
	}


	acl::VectorOfElements LevelSet::isBoundaryEl(unsigned int iEl)
	{
		unsigned int n(vto->elementaryCells[iEl].getSize());
		acl::VectorOfElements v(n);
		for(unsigned int i(0); i < n; ++i)
			v[i] = distanceTVE->values[vto->elementaryCells[iEl][i]];
		return acl::differentSign(v);
	}

	acl::VectorOfElements LevelSet::isBoundaryDir(unsigned int iDir)
	{
		return acl::differentSign(cat(subVE(distanceTVE->values,0), 
		                              subVE(distanceTVE->values,iDir))); 
	}

	acl::VectorOfElements LevelSet::getBoundaryCenter(unsigned int iEl)
	{
		auto cell(vto->elementaryCells[iEl]);
		auto &dVal(distanceTVE->values);
		acl::TypeID type(getElementType(dVal));
		unsigned int n(cell.getSize());
		vector<AVec<>> cellPoints;
		vto->getCellPoints(iEl, cellPoints);
		unsigned int nd(nD(cellPoints[0]));
		
		acl::VectorOfElements v(acl::generateVEConstantN(nd,0));
		acl::VectorOfElements nP(acl::generateVEConstant(0));
		for(unsigned int i(0); i < n; ++i)
			for(unsigned int j(i+1); j < n; ++j)
			{
				auto di(subVE(dVal, cell[i]));
				auto dj(subVE(dVal, cell[j]));

				auto isB(acl::differentSign(cat(di, dj)));
				auto x(di / (di - dj));
				auto pi(acl::generateVEConstant(cellPoints[i]));
				auto pj(acl::generateVEConstant(cellPoints[j]));

				auto p(select(acl::generateVEConstantN(nd,0), 
				              pi + (pj - pi) * x, 
				              isB,
				              type));
				copy(v + p, v);
				copy(nP +select(acl::generateVEConstant(0), 
				               acl::generateVEConstant(1), 
				               isB,
				               type),nP);
				
			}
		return v/nP;
	}

	vector<acl::Element> LevelSet::gcBoundaryArea(unsigned int iEl, 
	                                              acl::VectorOfElements & center,
	                                              acl::VectorOfElements & area)
	{
		auto cell(vto->elementaryCells[iEl]);
		auto &dVal(distanceTVE->values);
		acl::TypeID type(getElementType(dVal));
		unsigned int n(cell.getSize());
		vector<AVec<>> cellPoints;
		vto->getCellPoints(iEl, cellPoints);
//		unsigned int nd(nD(cellPoints[0]));

		vector<acl::Element> code(0);
		code.push_back((area = acl::generateVEConstant(0.))[0]);
		auto aij(acl::generateVEPrivateVariable(1u, type));
		auto aik(acl::generateVEPrivateVariable(1u, type));
		auto ajk(acl::generateVEPrivateVariable(1u, type));
		
		acl::VectorOfElements nP(acl::generateVEConstant(0));
		for(unsigned int i(0); i < n-2; ++i)
			for(unsigned int j(i+1); j < n-1; ++j)
				for(unsigned int k(j+1); k < n; ++k)
				{
					auto di(subVE(dVal, cell[i]));
					auto dj(subVE(dVal, cell[j]));
					auto dk(subVE(dVal, cell[k]));

					auto isBij(acl::differentSign(cat(di, dj)));
					auto isBik(acl::differentSign(cat(di, dk)));
					auto isBjk(acl::differentSign(cat(dj, dk)));
					auto xij(di / (di - dj));
					auto xik(di / (di - dk));
					auto xjk(dj / (dj - dk));
					auto pi(acl::generateVEConstant(cellPoints[i]));
					auto pj(acl::generateVEConstant(cellPoints[j]));
					auto pk(acl::generateVEConstant(cellPoints[k]));

					auto pij(pi + (pj - pi) * xij);
					auto pik(pi + (pk - pi) * xik);
					auto pjk(pj + (pk - pj) * xjk);
					auto nij(crossProduct(pik-center,pjk-center));
					auto nik(crossProduct(pij-center,pjk-center));
					auto njk(crossProduct(pij-center,pik-center));

					auto codeij(gcLength2(nij,aij));
					auto codeik(gcLength2(nij,aik));
					auto codejk(gcLength2(nij,ajk));
					code.reserve(code.size()+codeij.size()+codeik.size()+codejk.size());
					copy(codeij.begin(),codeij.end(),back_inserter(code));
					copy(codeik.begin(),codeik.end(),back_inserter(code));
					copy(codejk.begin(),codejk.end(),back_inserter(code));
					
					auto res(select(acl::generateVEConstant(0.),
					                   fabs(select(select(ajk, aik, !isBik), aij, !isBij)),
					                   isBij || isBik || isBjk, 
					                   type));
					code.push_back((area+=res)[0]);
				}
		return code;
	}
		
		
	acl::VectorOfElements LevelSet::getBoundaryPoint(unsigned int iDir)
	{
		auto &dVal(distanceTVE->values);
		
		auto d0(subVE(dVal, 0));
		auto di(subVE(dVal, iDir));

		auto x(d0 / (d0 - di));
		auto pi(acl::generateVEConstant(vectorTemplate->vectors[iDir]));
		return pi * x;
	}

		
	acl::VectorOfElements LevelSet::getValueOnBoundary(acl::VectorOfElements field,
	                                                   unsigned int i)
	{
		auto d0(subVE(distanceTVE->values, 0));
		auto di(subVE(distanceTVE->values, i));
		return subVE(field,0) + 
			(subVE(field,i) - subVE(field,0)) * d0 / (d0 - di);
	}
		
		
	void LevelSet::init()
	{
		distanceFieldInternalData=clone(distanceField);
		acl::initData(distanceFieldInternalData->getDContainer(),
		              distanceField->getEContainer());

		initKernelPropagation();		
	}

	void LevelSet::execute()
	{
		kernel->compute();
		swapBuffers(distanceField->getDContainer(),
		            distanceFieldInternalData->getDContainer());
	}
		
} // asl

