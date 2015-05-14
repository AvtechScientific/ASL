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


#include "aslBCond.h"
#include <acl/acl.h>
#include <acl/DataTypes/aclMemBlock.h>
#include <acl/aclGenerators.h>
//#include <agl/aglPointsList.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <math/aslTemplates.h>
#include <data/aslDataWithGhostNodes.h>
#include <aslGenerators.h>
#include <math/aslTemplateVE.h>
#include <math/aslDistanceFunctionAlg.h>

namespace asl
{

	BCond::BCond(const Block & b):
		block(b),
		templ(nearestNeigboursVT(nD(b)))
	{		
	}

	BCond::BCond(const Block & b, const VectorTemplate *const t):
		block(b),
		templ(t)
	{		
	}
		
	void BCond::addPoint(AVec<int> a,int d,double fr)
	{
		indices.push_back(block.c2i(a));
		directions.push_back(d);
		fractions.push_back(fr);
	}


	void BCond::loadIndicesToACL()
	{
		indicesACL=acl::SPVectorOfElementsData(new acl::VectorOfElementsData(1,indices.size(),int(0)));
		acl::copy(indices,(*indicesACL)[0]);
	}


	void BCond::loadDirectionsToACL()
	{
		directionsACL=acl::SPVectorOfElementsData(new acl::VectorOfElementsData(1,indices.size(),int(0)));
		acl::copy(directions,(*directionsACL)[0]);                       
	}


	void BCond::loadNeighbourIndicesToACL()
	{
		neighbourIndices.resize(indices.size());
		for (unsigned int i(0); i < indices.size(); ++i)
			neighbourIndices[i] =indices[i]+block.c2i(templ->vectors[directions[i]]);
			
		neighbourIndicesACL=acl::SPVectorOfElementsData(new acl::VectorOfElementsData(1,indices.size(),int(0)));
		acl::copy(neighbourIndices,(*neighbourIndicesACL)[0]);                       
	}

		
	const Block& BCond::getBlock()
	{
		return block;	
	}

	BCondWithMap::BCondWithMap(SPAbstractDataWithGhostNodes m, const VectorTemplate * const vt):
			pointsListFlag(false),
			currentPoint(acl::generateVEIndex()),
			templ(vt),
			bl(m->getBlock()),
			map(m)
	{
		unsigned int numD(nD(bl));
		if(numD != 2 && numD != 3)
			errorMessage("BCondWithMap : the map has wrong dimentionality");	

	}	

	BCondWithMap::BCondWithMap(SPDistanceFunction m, 
	                           const Block & b, 
	                           const VectorTemplate * const vt):
		pointsListFlag(false),
		currentPoint(acl::generateVEIndex()),
		templ(vt),
		bl(b),
		mapDF(m)
	{
		unsigned int numD(nD(bl));
		if(numD != 2 && numD != 3)
			errorMessage("BCondWithMap : the map has wrong dimentionality");	
	}	
		
	BCondWithMap::BCondWithMap(SPAbstractDataWithGhostNodes m, 
	                           SPAbstractDataWithGhostNodes cd, 
	                           const VectorTemplate * const vt):
		pointsListFlag(false),
		currentPoint(acl::generateVEIndex()),
		templ(vt),
		bl(m->getBlock()),
		map(m),
		computationalDomain(cd)	
	{
		unsigned int numD(nD(bl));
		unsigned int numDcd(nD(computationalDomain->getBlock()));
		if(numD != numDcd)
			errorMessage("BCondWithMap : the map and computationalDomain has different dimentionality");	
		if(numD != 2 && numD != 3)
			errorMessage("BCondWithMap : the map has wrong dimentionality");	

	}				

	BCondWithMap::BCondWithMap(SPAbstractDataWithGhostNodes m, 
	                           SPDistanceFunction cd, 
	                           const VectorTemplate * const vt):
		pointsListFlag(false),
		currentPoint(acl::generateVEIndex()),
		templ(vt),
		bl(m->getBlock()),
		map(m),
		computationalDomainDF(cd)	
	{
		unsigned int numD(nD(bl));
		if(numD != 2 && numD != 3)
			errorMessage("BCondWithMap : the map has wrong dimentionality");	
	}				

	BCondWithMap::BCondWithMap(SPDistanceFunction m, 
	                           SPDistanceFunction cd,
	                           const Block & b,
	                           const VectorTemplate * const vt):
		pointsListFlag(false),
		currentPoint(acl::generateVEIndex()),
		templ(vt),
		bl(b),
		mapDF(m),
		computationalDomainDF(cd)	
	{
		unsigned int numD(nD(bl));
		if(numD != 2 && numD != 3)
			errorMessage("BCondWithMap : the map has wrong dimentionality");	
	}				
		
		
	void BCondWithMap::initMapInfrastructure(acl::ExpressionContainer & ec)
	{
		bool isMapDF(mapDF.get()!=0);
		bool isCDDF(computationalDomainDF.get()!=0);
		//todo!!!!!!!!!!!!!!!!!!
		if(isMapDF || isCDDF)
		{
			
		}
			
		if(isMapDF)
		{
//			mapTVE.reset(new TemplateVE(*mapDF,, *templ));
			ec << mapTVE->initValues;		
			
		}
		else
		{
			auto mapX(generateDCFullSafe(map, -.9));
			mapTVE.reset(new TemplateVE(*mapX, *templ));
//			mapTVE.reset(new TemplateVE(*map, *templ));
			ec << mapTVE->initValues;		
		}
		bool isCompDom(computationalDomain.get()!=0);
		if(isCompDom)
		{
			auto compDomX(generateDCFullSafe(computationalDomain, -1));
			cDomainTVE.reset(new TemplateVE(*compDomX, *templ));
			ec << cDomainTVE->initValues;
		}
	}

	acl::VectorOfElements BCondWithMap::isGhostNode(unsigned int i)
	{
		return asl::isGhostNode(*mapTVE,i);
	}

	acl::VectorOfElements BCondWithMap::isComputationNode(unsigned int i)
	{
		auto res(asl::isComputationNode(*mapTVE,i));
		bool isCompDom(computationalDomain.get()!=0);
		if(isCompDom)
			copy(res && asl::isComputationNode(*cDomainTVE,i), res);
		return res;
	}

	acl::VectorOfElements BCondWithMap::isComputationNode(const vector<unsigned int> & ii)
	{
		auto res(subVE(mapTVE->values,ii[0]) > 0.);
		for(unsigned int i(1); i < ii.size(); ++i)
			copy(res && (subVE(mapTVE->values,ii[i]) > 0.), res);
		bool isCompDom(computationalDomain.get()!=0);
		if(isCompDom)
		for(unsigned int i(0); i < ii.size(); ++i)			
			copy(res && (acl::subVE(cDomainTVE->values,ii[i]) > 0.), res);
		return res;
	}
		
	acl::VectorOfElements BCondWithMap::isGhostNode()
	{
		return (map->getEContainer() <= 0.);
	}

	acl::VectorOfElements BCondWithMap::isComputationNode()
	{
		auto res(map->getEContainer() > 0.);
		bool isCompDom(computationalDomain.get()!=0);
		if(isCompDom)
			copy(res && (computationalDomain->getEContainer() > 0.), res);
		return res;
	}
		
	BCondSlice::BCondSlice(const Block & b):
		block(b),
		templ(NULL)
	{		
	}

	BCondSlice::BCondSlice(const Block & b, const VectorTemplate *const t):
		block(b),
		templ(t)
	{		
	}

	inline unsigned int count0(const AVec<int> &a)
	{
		unsigned int n(0);
		for (unsigned int i(0); i < nD(a); ++i)
			n+=a[i] ==0?1:0;
		return n;
	}
		
	void BCondSlice::addGhostSlice(AVec<int> pB,AVec<int> pE,int d)
	{
		pointB=block.c2i(pB);
		pointE=block.c2i(pE);
		direction=d;
		AVec<int> dif (pE-pB);
		unsigned int n(count0(dif));
		sliceDimentions.resize(n);
		sliceIncrements.resize(n);
		for (unsigned int i(0),k(0); i < nD(dif); ++i){
			sliceDimentions[k] =dif[i];
			sliceIncrements[k] =block.c2iTransformVector[i];
			k+=dif[i] ==0?0:1;
		}
	}

	const Block& BCondSlice::getBlock()
	{
		return block;	
	}
		

	BCondConnector::BCondConnector(const Block & b1, const Block & b2):
		block1(b1),
		block2(b2)
	{		
	}


	void BCondConnector::addGhostPoint(AVec<int> a1,AVec<int> a2)
	{
		indices1.push_back(block1.c2i(a1));
		directions1.push_back(0);
		indices2.push_back(block2.c2i(a2));
		directions2.push_back(0);
	}


	void BCondConnector::addGhostPoint(AVec<int> a1,int d1,AVec<int> a2,int d2)
	{
		indices1.push_back(block1.c2i(a1));
		directions1.push_back(d1);
		indices2.push_back(block2.c2i(a2));
		directions2.push_back(d2);
	}

	const Block & BCondConnector::getBlock1()
	{
		return block1;	
	}

	const Block & BCondConnector::getBlock2()
	{
		return block2;	
	}


	void BCondConnector::loadIndicesToACL()
	{
		indices1ACL=acl::SPVectorOfElementsData(new acl::VectorOfElementsData(1,indices1.size(),int(0)));
		acl::copy(indices1,(*indices1ACL)[0]);                       
		indices2ACL=acl::SPVectorOfElementsData(new acl::VectorOfElementsData(1,indices2.size(),int(0)));
		acl::copy(indices2,(*indices2ACL)[0]);                       
	}


	void BCondConnector::loadDirectionsToACL()
	{
		directions1ACL=acl::SPVectorOfElementsData(new acl::VectorOfElementsData(1,indices1.size(),int(0)));
		acl::copy(directions1,(*directions1ACL)[0]);                       
		directions2ACL=acl::SPVectorOfElementsData(new acl::VectorOfElementsData(1,indices2.size(),int(0)));
		acl::copy(directions2,(*directions2ACL)[0]);                       
	}

	BCondConnectorSlice::BCondConnectorSlice(const Block & b1, const Block & b2, const VectorTemplate *const t):
		block1(b1),
		block2(b2),
		templ(t)
	{		
	}
		
	void BCondConnectorSlice::addGhostSlice1(AVec<int> pB,AVec<int> pE,int d)
	{
		point1B=block1.c2i(pB);
		point1E=block1.c2i(pE);
		direction1=d;
		AVec<int> dif (pE-pB);
		unsigned int n(count0(dif));
		sliceDimentions1.resize(n);
		sliceIncrements1.resize(n);
		for (unsigned int i(0),k(0); i < nD(dif); ++i){
			sliceDimentions1[k] =dif[i];
			sliceIncrements1[k] =block1.c2iTransformVector[i];
			k+=dif[i] ==0?0:1;
		}
	}


	void BCondConnectorSlice::addGhostSlice2(AVec<int> pB,AVec<int> pE,int d)
	{
		point2B=block1.c2i(pB);
		point2E=block1.c2i(pE);
		direction2=d;
		AVec<int> dif (pE-pB);
		unsigned int n(count0(dif));
		sliceDimentions2.resize(n);
		sliceIncrements2.resize(n);
		for (unsigned int i(0),k(0); i < nD(dif); ++i){
			sliceDimentions2[k] =dif[i];
			sliceIncrements2[k] =block2.c2iTransformVector[i];
			k+=dif[i] ==0?0:1;
		}
	}
		
	void addSliceX0(BCond & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		const VectorTemplate * vt(a.getVT());
		if (nD(s)==2){
			AVec<int> dir(vt->vectors[1]);
			for (unsigned int k(1); k<vt->vectors.size(); ++k)
				if(vt->vectors[k]*dir>0)
					for (int i(1); i < s[1]-1; ++i)
						a.addPoint(makeAVec (1,i)-vt->vectors[k],k);
		}
		else
			if (nD(s)==3){
				AVec<int> dir(vt->vectors[1]);
				for (unsigned int k(1); k<vt->vectors.size(); ++k)
					if(vt->vectors[k]*dir>0)
						for (int i(1); i < s[1]-1; ++i)
							for (int j(1);j<s[2]-1;++j)
								a.addPoint(makeAVec (1,i,j)-vt->vectors[k],k);
			}			
	}


	void addSliceXE(BCond & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		const VectorTemplate * vt(a.getVT());
		if (nD(s)==2){
			AVec<int> dir(vt->vectors[3]);
			for (unsigned int k(1); k<vt->vectors.size(); ++k)
				if(vt->vectors[k]*dir>0)
					for (int i(1); i < s[1]-1; ++i)
						a.addPoint(makeAVec (s[0]-2,i)-vt->vectors[k],k);
		}
		else
			if (nD(s)==3){
				AVec<int> dir(vt->vectors[4]);
				for (unsigned int k(1); k<vt->vectors.size(); ++k)
				if(vt->vectors[k]*dir>0)
					for (int i(1); i < s[1]-1; ++i)
						for (int j(1);j<s[2]-1;++j)
							a.addPoint(makeAVec (s[0]-2,i,j)-vt->vectors[k],k);
			}			
	}


	void addSliceY0(BCond & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		const VectorTemplate * vt(a.getVT());
		if (nD(s)==2){
			AVec<int> dir(vt->vectors[2]);
			for (unsigned int k(1); k<vt->vectors.size(); ++k)
				if(vt->vectors[k]*dir>0)
					for (int i(1); i < s[0]-1; ++i)
						a.addPoint(makeAVec (i,1)-vt->vectors[k],k);
		}
		else
			if (nD(s)==3){
				AVec<int> dir(vt->vectors[2]);
				for (unsigned int k(1); k<vt->vectors.size(); ++k)
					if(vt->vectors[k]*dir>0)
						for (int i(1); i < s[0]-1; ++i)
							for (int j(1);j < s[2]-1;++j)
								a.addPoint(makeAVec (i,1,j)-vt->vectors[k],k);
			}			
			
	}


	void addSliceYE(BCond & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		const VectorTemplate * vt(a.getVT());
		if (nD(s)==2){
			AVec<int> dir(vt->vectors[4]);
			for (unsigned int k(1); k<vt->vectors.size(); ++k)
				if(vt->vectors[k]*dir>0)
					for (int i(1); i < s[0]-1; ++i)
						a.addPoint(makeAVec (i,s[1]-2)-vt->vectors[k],k);
		}else
			if (nD(s)==3){
				AVec<int> dir(vt->vectors[5]);
				for (unsigned int k(1); k<vt->vectors.size(); ++k)
					if(vt->vectors[k]*dir>0)
						for (int i(1); i < s[0]-1; ++i)
							for (int j(1); j < s[2]-1; ++j)
								a.addPoint(makeAVec (i,s[1]-2,j)-vt->vectors[k],k);
			}			
	}


	void addSliceZ0(BCond & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		const VectorTemplate * vt(a.getVT());
		if (nD(s)==3){
			AVec<int> dir(vt->vectors[3]);
			for (unsigned int k(1); k<vt->vectors.size(); ++k)
				if(vt->vectors[k]*dir>0)
					for (int i(1); i < s[0]-1; ++i)
						for (int j(1); j < s[1]-1; ++j)
							a.addPoint(makeAVec (i,j,1)-vt->vectors[k],k);
		}			
			
	}


	void addSliceZE(BCond & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		const VectorTemplate * vt(a.getVT());
		if (nD(s)==3){
			AVec<int> dir(vt->vectors[6]);
			for (unsigned int k(1); k<vt->vectors.size(); ++k)
				if(vt->vectors[k]*dir>0)
					for (int i(1); i < s[0]-1; ++i)
						for (int j(1);j<s[1]-1;++j)
							a.addPoint(makeAVec (i,j,s[2]-2)-vt->vectors[k],k);
		}			
	}

	void addSlices(BCond & bc, const vector<SlicesNames> & sl)
	{
		if(find(sl.begin(), sl.end(), X0) != sl.end())
			addSliceX0(bc);
		if(find(sl.begin(), sl.end(), XE) != sl.end())
			addSliceXE(bc);
		if(find(sl.begin(), sl.end(), Y0) != sl.end())
			addSliceY0(bc);
		if(find(sl.begin(), sl.end(), YE) != sl.end())
			addSliceYE(bc);
		if(find(sl.begin(), sl.end(), Z0) != sl.end())
			addSliceZ0(bc);
		if(find(sl.begin(), sl.end(), ZE) != sl.end())
			addSliceZE(bc);
		
	}

		
	void addSliceX(BCond & a, int x)
	{
		const AVec<int> & s(a.getBlock().getSize());
		if (nD(s)==2){
			for (int i(0); i < s[1]; ++i)
				a.addPoint(makeAVec (x,i));
		}
		else
			if (nD(s)==3){
				for (int i(0); i < s[1]; ++i)
					for (int j(0);j<s[2];++j)
						a.addPoint(makeAVec (x,i,j));
			}		
	}

	void addSliceY(BCond & a, int y)
	{
		const AVec<int> & s(a.getBlock().getSize());
		if (nD(s)==2){
			for (int i(0); i < s[0]; ++i)
				a.addPoint(makeAVec (i,y));
		}
		else
			if (nD(s)==3){
				for (int i(0); i < s[0]; ++i)
					for (int j(0);j<s[2];++j)
						a.addPoint(makeAVec (i,y,j));
			}		
	}

	void addSliceZ(BCond & a, int z)
	{
		const AVec<int> & s(a.getBlock().getSize());
		if (nD(s)==2){
			asl::errorMessage("addSliceZ: The data dimentionality is 2");
		}
		else
			if (nD(s)==3){
				for (int i(0); i < s[0]; ++i)
					for (int j(0);j<s[1];++j)
						a.addPoint(makeAVec (i,j,z));
			}		
	}
		

	void addSliceX0(BCondSlice & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		if (nD(s)==2)
			a.addGhostSlice(makeAVec(0,0),makeAVec(0,s[1]),1);
		else
			if (nD(s)==3)
				a.addGhostSlice(makeAVec(0,0,0),makeAVec(0,s[1],s[2]),1);
	}


	void addSliceXE(BCondSlice & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		if (nD(s)==2)
			a.addGhostSlice(makeAVec(s[0],0),makeAVec(s[0],s[1]),3);
		else
			if (nD(s)==3)
				a.addGhostSlice(makeAVec(s[0],0,0),makeAVec(s[0],s[1],s[2]),4);
	}


	void addSliceY0(BCondSlice & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		if (nD(s)==2)
			a.addGhostSlice(makeAVec(0,0),makeAVec(s[0],0),2);
		else
			if (nD(s)==3)
				a.addGhostSlice(makeAVec(0,0,0),makeAVec(s[0],0,s[2]),2);
	}


	void addSliceYE(BCondSlice & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		if (nD(s)==2)
			a.addGhostSlice(makeAVec(0,s[1]),makeAVec(s[0],s[1]),4);
		else
			if (nD(s)==3)
				a.addGhostSlice(makeAVec(0,s[1],0),makeAVec(s[0],s[1],s[2]),5);
	}


	void addSliceZ0(BCondSlice & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		if (nD(s)<3)
			errorMessage ("addSliceZ0: The block dimensionality is less than 3");
		else
			if (nD(s)==3)
				a.addGhostSlice(makeAVec(0,0,0),makeAVec(s[0],s[1],0),3);
	}


	void addSliceZE(BCondSlice & a)
	{
		const AVec<int> & s(a.getBlock().getSize());
		if (nD(s)<3)
			errorMessage ("addSliceZE: The block dimensionality is less than 3");
		else
			if (nD(s)==3)
				a.addGhostSlice(makeAVec(0,0,s[2]),makeAVec(s[0],s[1],s[2]),6);
	}

	/// checks intersection of a litice vector with a plane
	/**
		\param p0 point on the plane
		\param n  plane normal
		\param pn current node
		\param templ the used template
		\param dir number of direction
	*/
	inline bool checkIntersection(const AVec<> & p0, 
	                               const AVec<> & n, 
	                               const AVec<> & pn,
	                               const VectorTemplate * templ, 
	                               unsigned int dir)
	{
		double x(((p0-pn)*n)/(n*AVec<>(templ->vectors[dir])));
		return x>0 && x<1; 
	}

	/// adds a point to a BCond object
	/**
		\param p0 point on the plane
		\param n  plane normal
		\param pn current node
		\param templ the used template
		\param dir number of direction
	*/
	inline void addPointToBC(BCond  & bc, 
	             			  const AVec<> & p0, 
	             			  const AVec<> & n, 
	             			  const AVec<int> & pn)
	{
		if (in(bc.getBlock(),pn))
		{
			AVec<> r(AVec<>(pn)*bc.getBlock().dx+bc.getBlock().position);
			double dist((r-p0)*n);
			if (dist<=0 && dist>2.)
				for(unsigned int i(0); i<bc.getVT()->vectors.size();++i)
					if (checkIntersection(p0,n,r,bc.getVT(),i))
						bc.addPoint(pn,i);
		}
		
	}
		
/*	void addAGLObject(BCond  & bc, agl::SPPointsList geom)
	{
		unsigned int maxId(geom->maxId());		
		for (unsigned int i(0); i < maxId; geom->next(i))
		{
			AVec<> p0(geom->getPoint(i));
			AVec<int> pn0(asl::round((p0-bc.getBlock().position)/bc.getBlock().dx));
			for(unsigned int i(0); i<bc.getVT()->vectors.size(); ++i)
			{
				AVec<int> pn(pn0+bc.getVT()->vectors[i]);
				addPointToBC(bc, geom->getPoint(i), geom->getSurfaceNormal(i), pn);
			}
		}
	}
*/		
} // asl

