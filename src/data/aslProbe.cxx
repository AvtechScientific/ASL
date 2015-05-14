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


#include "aslProbe.h"
#include "aslDataWrapper.h"
#include <acl/acl.h>
#include <acl/aclGenerators.h>
#include <acl/DataTypes/aclArray.h>
#include <acl/aclMath/aclVectorOfElements.h>
#include <math/aslTemplateVE.h>
#include <aslGenerators.h>


namespace asl
{

	Probe::Probe(SPAbstractData d):
		data(d),
		values(getNComponents())
	{		
	}

	void Probe::addPoint(AVec<int> p)
	{
		if (p.getSize() != nD(data->getBlock()))
			errorMessage("Probe::addPoint() - attempt to add point that has wrong dimensions");

		if (!in(data->getBlock(), p))
			errorMessage("Probe::addPoint() - attempt to add point that is out of block range");
		
		indices.push_back(data->getBlock().c2i(p));
	}


	const unsigned int Probe::getNComponents() const
	{
		return data->getEContainer().size();		
	}


	const unsigned int Probe::getNDimensions() const
	{
		return nD(data->getBlock());
	}

		
	void Probe::loadIndicesToACL()
	{
		indicesACL = acl::SPVectorOfElementsData(new acl::VectorOfElementsData(1,
		                                                                       indices.size(),
		                                                                       int(0)));
		acl::copy(indices, (*indicesACL)[0]);
	}


	void Probe::loadValuesFromACL()
	{
		unsigned int nC(getNComponents());
		for (unsigned int i(0); i < nC; ++i)
			acl::copy((*valuesACL)[i], values[i]);
	}
	

	void Probe::init()
	{
		loadIndicesToACL();
		valuesACL = acl::SPVectorOfElementsData(new acl::VectorOfElementsData(getNComponents(),
		                                                                      indices.size(),
		                                                                      double(0)));	
		k << (*valuesACL = excerpt(data->getEContainer(), *indicesACL));
		k.setup();
	}


	void Probe::update()
	{
		k.compute();
		loadValuesFromACL();
	}


	ProbeLI::ProbeLI(SPAbstractData d):
		data(d),
		values(getNComponents())
	{		
	}
		
	void ProbeLI::addPoint(AVec<> p)
	{
		if (p.getSize() != nD(data->getBlock()))
			errorMessage("Probe::addPoint() - attempt to add point that has wrong dimensions");

		if (!in(data->getBlock(), p))
			errorMessage("Probe::addPoint() - attempt to add point that is out of block range");
		
		points.push_back(p);
	}


	const unsigned int ProbeLI::getNComponents() const
	{
		return data->getEContainer().size();		
	}


	const unsigned int ProbeLI::getNDimensions() const
	{
		return nD(data->getBlock());
	}
		
	void ProbeLI::loadPointsToACL()
	{
		unsigned int nd(getNDimensions());
		unsigned int length(points.size());
		pointsACL = acl::SPVectorOfElementsData(new acl::VectorOfElementsData(nd, length, double(0)));
		vector<vector<double>> pointsCompon(nd);
		for(unsigned int i(0); i<nd; ++i)
			pointsCompon[i].resize(points.size());

		for(unsigned int j(0); j<length; ++j)
			for(unsigned int i(0); i<nd; ++i)
				pointsCompon[i][j]=points[j][i];

		for(unsigned int i(0); i<nd; ++i)
			acl::copy(pointsCompon[i], (*pointsACL)[i]);
	}


	void ProbeLI::loadValuesFromACL()
	{
		unsigned int nC(getNComponents());
		for (unsigned int i(0); i < nC; ++i)
			acl::copy((*valuesACL)[i], values[i]);
	}
	

	void ProbeLI::init()
	{
		loadPointsToACL();
		auto & bl(data->getBlock());
		unsigned int nd(pointsACL->size()); 
		valuesACL = acl::SPVectorOfElementsData(new acl::VectorOfElementsData(getNComponents(),
		                                                                      points.size(),
		                                                                      double(0)));
		auto ind(acl::generateVEPrivateVariable<int>(1));
		auto e(acl::generateVEPrivateVariable<double>(nd));

		k << (e = (*pointsACL - bl.position)/bl.dx);
		k << (ind = floor(e) * bl.c2iTransformVector);
		k << (e = e - acl::floor(e));

		TemplateVE dTVE;
		auto d(generateDataContainer_SP(data->getBlock(), data->getEContainer(),0));
		for(unsigned int i(0); i < valuesACL->size(); ++i)
		{
			dTVE.init(*d, *elementaryCellVT(nd), i);
			k << excerpt(dTVE.initValues,ind);
			k << (subVE(*valuesACL,i) = interpolate(dTVE, e));
		}
		k.setup();
	}


	void ProbeLI::update()
	{
		k.compute();
		loadValuesFromACL();
	}

		
}// namespace asl
