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


#ifndef ASLPROBETEMPLATES_H
#define ASLPROBETEMPLATES_H

#include "aslVectors.h"

namespace asl {


	class ProbeTemplate
	{
		public:
			std::vector<AVec<double> > vectors;
			ProbeTemplate(int n, AVec<double>* vec);		

			inline unsigned int numberOfDimentions() const;
	};

	
	/// A triangular probe 1D space
	/**	
	 \ingroup Templates
		 
	 \image html probeD1Q2.png "1D triangular probe"
	*/
	extern ProbeTemplate probeD1q2;

	
	/// A triangular probe 2D space
	/**	
	 \ingroup Templates
		 
	 \image html probeD2Q3.png "2D triangular probe"
	*/
	extern ProbeTemplate probeD2q3;

	/// A triangular probe 3D space
	/**	
	 \ingroup Templates
		 
	 \image html probeD3Q4.png "3D triangular probe"
	*/
	extern ProbeTemplate probeD3q4;


	/// returns template corresponding to minimal probes \ingroup Templates
	/** asl::probeD1q2, asl::probeD2q3, asl::probeD3q4
	*/
	inline const ProbeTemplate* allMinimalProbeTemplate(unsigned int dimNumber);
	
// ----------------------------- Implementation -------------------------

	unsigned int ProbeTemplate::numberOfDimentions() const
	{
		return vectors[0].getSize();
	}

	inline const ProbeTemplate* minimalProbeTemplate(unsigned int dimNumber)
	{
		static ProbeTemplate* vt[3]={&probeD1q2,&probeD2q3,&probeD3q4};
		return vt[dimNumber-1];
	}
	
}// asl

#endif // TEMPL_H_INCLUDED
