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


#include "aslDistanceFunction.h"
#include "acl/aclMath/aclVectorOfElements.h"
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include <data/aslBlocks.h>
#include <aslGenerators.h>
#include <math/aslIndex2Position.h>
#include "acl/aclMath/aclMathAlg.h"
#include <math/aslTemplatesExtras.h>
#include <math/aslTemplateVE.h>



namespace asl
{

	acl::VectorOfElements isGhostNode(TemplateVE & distanceTVE, unsigned int i)
	{
		return (distanceTVE.getValue(i) <= 0.);
	}

	acl::VectorOfElements isComputationNode(TemplateVE & distanceTVE, unsigned int i)
	{
		return (distanceTVE.getValue(i) > 0.);
	}


	acl::VectorOfElements nGhostNodesInCell(TemplateVE & distanceTVE, unsigned int iEl)
	{
		auto cell(distanceTVE.vto->elementaryCells[iEl]);
		acl::VectorOfElements nGhost(1);
		acl::TypeID type(getElementType(distanceTVE.values));
		copy(select(acl::generateVEConstant(1), 
		            isGhostNode(distanceTVE,cell[0]),
		            type),
		     nGhost);
			
		for(unsigned int i(1); i< cell.getSize(); ++i)
		{
			copy(nGhost + select(acl::generateVEConstant(1), 
			                     isGhostNode(distanceTVE,cell[i]),
			                     type),
			     nGhost);
		}
		return nGhost;
	}

	acl::VectorOfElements isBoundaryBetween(TemplateVE & distanceTVE, 
	                                        unsigned int iE, 
	                                        unsigned int i, 
	                                        unsigned int j)
	{
		auto vto(distanceTVE.vto);
		auto cell(vto->elementaryCells[iE]);
		return (isGhostNode(distanceTVE, cell[i]) && 
		        isComputationNode(distanceTVE, cell[j])) || 
			   (isGhostNode(distanceTVE, cell[j]) && 
			    isComputationNode(distanceTVE, cell[i])); 		
	}
	
	acl::VectorOfElements exBoundaryX(TemplateVE & distanceTVE, unsigned int i)
	{
		return distanceTVE.getValue(0)/(distanceTVE.getValue(0)-distanceTVE.getValue(i));
	}
	
	acl::VectorOfElements exBoundaryCenter(TemplateVE & distanceTVE, unsigned int iEl)
	{
		auto vto(distanceTVE.vto);
		auto cell(vto->elementaryCells[iEl]);
		auto &dVal(distanceTVE.values);
		acl::TypeID type(getElementType(dVal));
		unsigned int n(cell.getSize());
		vector<AVec<>> cellPoints;
		vto->getCellPoints(iEl, cellPoints);
		unsigned int nd(nD(cellPoints[0]));
		
		acl::VectorOfElements v(acl::generateVEConstantN(nd,0));
		acl::VectorOfElements nP(acl::generateVEConstant(0));
		for(unsigned int i(0); i < n-1; ++i)
			for(unsigned int j(i+1); j < n; ++j)
			{
				auto di(subVE(dVal, cell[i]));
				auto dj(subVE(dVal, cell[j]));

				auto isB(isBoundaryBetween(distanceTVE,iEl,i,j));
				auto x(di / (di - dj));
				auto pi(acl::generateVEConstant(cellPoints[i]));
				auto pj(acl::generateVEConstant(cellPoints[j]));

				auto p(select(pi + (pj - pi) * x, isB, type));
				copy(v + p, v);
				copy(nP +select(acl::generateVEConstant(1), isB, type),nP);
				
			}
		return v/nP;
	}

	vector<acl::Element> gcBoundaryArea(TemplateVE & distanceTVE,
	                                    unsigned int iEl, 
	                                    acl::VectorOfElements & center,
	                                    acl::VectorOfElements & area)
	{
		auto vto(distanceTVE.vto);

		auto cell(vto->elementaryCells[iEl]);
		auto &dVal(distanceTVE.values);
		acl::TypeID type(getElementType(area));
		unsigned int n(cell.getSize());
		vector<AVec<>> cellPoints;
		vto->getCellPoints(iEl, cellPoints);
//		unsigned int nd(nD(cellPoints[0]));

		vector<acl::Element> code(0);
		code<<(area = acl::generateVEConstant(0.));
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

					auto isBij(isBoundaryBetween(distanceTVE,iEl, i, j));
					auto isBik(isBoundaryBetween(distanceTVE,iEl, i, k));
					auto isBjk(isBoundaryBetween(distanceTVE,iEl, j, k));
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
					auto codeik(gcLength2(nik,aik));
					auto codejk(gcLength2(njk,ajk));
					code<<codeij;
					code<<codeik;
					code<<codejk;

					auto res(select(sqrt(select(select(ajk, aik, !isBik), aij, !isBij)),
					                isBij || isBik || isBjk, type)*.5);
					code<<(area+=res);
				}
		return code;
	}

	vector<acl::Element> gcBoundaryArea(TemplateVE & distanceTVE,
	                                    acl::VectorOfElements & area)
	{
		unsigned int nCells(distanceTVE.vto->elementaryCells.size());
		acl::TypeID type(getElementType(area));
		
		vector<acl::Element> code;
		auto areaCell(acl::generateVEPrivateVariable(1u, type));
		auto centerCell(acl::generateVEPrivateVariable(3u, type));

		code<<(area = acl::generateVEConstant(0.));
		
		for(unsigned int i(0); i<nCells; ++i)
		{
			code<<(centerCell = exBoundaryCenter(distanceTVE, i));
			code<<gcBoundaryArea(distanceTVE,i,centerCell,areaCell);
			code<<(area+=areaCell);
		}
		

		return code;
	}

	acl::VectorOfElements zeroNodeWeight(TemplateVE & distanceTVE, unsigned int iEl)
	{
		auto cell(distanceTVE.vto->elementaryCells[iEl]);
		acl::VectorOfElements w0(1);
		acl::VectorOfElements wC(1);
		acl::TypeID type(getElementType(distanceTVE.values));
		copy(select(acl::generateVEConstant(distanceTVE.vectorTemplate->laplasCoefs[cell[1]]), 
		            isComputationNode(distanceTVE,cell[1]),
		            type),
		     w0);
			
		for(unsigned int i(2); i< cell.getSize(); ++i)
		{
			copy(w0+select(acl::generateVEConstant(distanceTVE.vectorTemplate->laplasCoefs[cell[i]]), 
				        isComputationNode(distanceTVE,cell[i]),
				        type),
				 w0);
		}
		return w0;
	}

	acl::VectorOfElements surfaceCellWeight(TemplateVE & distanceTVE, unsigned int iEl)
	{
		auto cell(distanceTVE.vto->elementaryCells[iEl]);
		acl::VectorOfElements w0(1);
		acl::VectorOfElements wC(1);
		acl::TypeID type(getElementType(distanceTVE.values));
		copy(select(acl::generateVEConstant(distanceTVE.vectorTemplate->laplasCoefs[cell[1]]), 
		            isComputationNode(distanceTVE,cell[1]),
		            type),
		     w0);
			
		for(unsigned int i(2); i< cell.getSize(); ++i)
		{
			copy(w0+select(acl::generateVEConstant(distanceTVE.vectorTemplate->laplasCoefs[cell[i]]), 
				        isComputationNode(distanceTVE,cell[i]),
				        type),
				 w0);
		}
		copy(w0,wC);
		for(unsigned int i(1); i<cell.getSize()-1; ++i)
			for(unsigned int j(i+1); j< cell.getSize(); ++j)
		{
			copy(wC+select(acl::generateVEConstant(edgeWeight(*distanceTVE.vto,iEl,i,j)), 
				        isBoundaryBetween(distanceTVE,iEl,i,j),
				        type),
				 wC);
		}
		
		return wC;
	}	
	vector<acl::Element> gcBoundaryAreaPerGhostPoint(TemplateVE & distanceTVE,
	                                                 unsigned int iEl, 
	                                                 acl::VectorOfElements & center, 
	                                                 acl::VectorOfElements & area)
	{
		auto code(gcBoundaryArea(distanceTVE,iEl,center,area));
		acl::TypeID type(getElementType(area));
//		auto wc(acl::generateVEPrivateVariable(1u, type));
		auto nGhost(acl::generateVEPrivateVariable(1u, type));
		code << (nGhost = nGhostNodesInCell(distanceTVE, iEl));
		
//		code << (wc = surfaceCellWeight(distanceTVE, iEl));
//		code << (area *=select(acl::generateVEConstant(1.), 
//		                       zeroNodeWeight(distanceTVE,iEl)/wc,
//		                       wc>.01,
//		                       type));
/*		code << (area *=select(acl::generateVEConstant(1.), 
		                       zeroNodeWeight(distanceTVE,iEl),
		                       wc>.01,
		                       type));
*/
		code << (area /=select(acl::generateVEConstant(1.), 
		                       nGhost,
		                       nGhost>.1,
		                       type));

		return code;
	}

	vector<acl::Element> gcBoundaryAreaPerGhostPoint(TemplateVE & distanceTVE,
	                                                 acl::VectorOfElements & area)
	{
		unsigned int nCells(distanceTVE.vto->elementaryCells.size());
		acl::TypeID type(getElementType(area));
		
		vector<acl::Element> code;
		auto areaCell(acl::generateVEPrivateVariable(1u, type));
		auto centerCell(acl::generateVEPrivateVariable(3u, type));

		code<<(area = acl::generateVEConstant(0.));
		
		for(unsigned int i(0); i<nCells; ++i)
		{
			code<<(centerCell = exBoundaryCenter(distanceTVE, i));
			code<<gcBoundaryAreaPerGhostPoint(distanceTVE,i,centerCell,areaCell);
			code<<(area+=areaCell);
		}

		return code;
	}
	
} // asl
