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


#ifndef ASLPROBE_H
#define ASLPROBE_H

#include <acl/Kernels/aclKernel.h>
#include <math/aslVectors.h>
#include <iostream>
#include <fstream>



using  namespace std;

namespace acl
{
	class VectorOfElementsData; 
	typedef std::shared_ptr<VectorOfElementsData> SPVectorOfElementsData;

}
namespace asl
{
	
	class AbstractData;
	typedef std::shared_ptr<AbstractData> SPAbstractData;
	
	/// Collects point values from the input data 
	/**
		/// \ingroup DataAnalysis
		 The class takes in
	*/
	class Probe
	{
		protected:
			SPAbstractData data;
			acl::Kernel k;

			std::vector<int> indices;
			acl::SPVectorOfElementsData indicesACL;

			vector<vector<double>> values;
			acl::SPVectorOfElementsData valuesACL;

			void loadIndicesToACL();
			void loadValuesFromACL();

		public:	
			Probe(SPAbstractData d);			
			void addPoint(AVec<int> p);
			/// initialization of internal kernels and data
			void init();
			/// Updates data in output
			void update();
			const unsigned int getNComponents() const;
			const unsigned int getNDimensions() const;
			inline vector<double> & getComponent(unsigned int i);
			inline AVec<double> getValue(unsigned int pointNumber);
	};	

	/// Collects point values from the input data with linear interpolation  
	/**
		/// \ingroup DataAnalysis
	*/
	class ProbeLI
	{
		protected:
			SPAbstractData data;
			acl::Kernel k;

			std::vector<AVec<>> points;
			acl::SPVectorOfElementsData pointsACL;

			vector<vector<double>> values;
			acl::SPVectorOfElementsData valuesACL;

			void loadPointsToACL();
			void loadValuesFromACL();

		public:	
			ProbeLI(SPAbstractData d);			
			void addPoint(AVec<> p);
			/// initialization of internal kernels and data
			void init();
			/// Updates data in output
			void update();
			const unsigned int getNComponents() const;
			const unsigned int getNDimensions() const;
			inline vector<double> & getComponent(unsigned int i);
			inline AVec<double> getValue(unsigned int pointNumber);
	};	


	// -------------------------- Implementation ------------------------------

	inline vector<double> & Probe::getComponent(unsigned int i)
	{
		return values[i];
	}


	inline AVec<double> Probe::getValue(unsigned int pointNumber)
	{
		unsigned int nC(getNComponents());
		AVec<double> value(nC);
		for (unsigned int i(0); i < nC; ++i)
			value[i] = values[i][pointNumber];

		return value;
	}

	inline vector<double> & ProbeLI::getComponent(unsigned int i)
	{
		return values[i];
	}


	inline AVec<double> ProbeLI::getValue(unsigned int pointNumber)
	{
		unsigned int nC(getNComponents());
		AVec<double> value(nC);
		for (unsigned int i(0); i < nC; ++i)
			value[i] = values[i][pointNumber];

		return value;
	}
	
}

#endif

