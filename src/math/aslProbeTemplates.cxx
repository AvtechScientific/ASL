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


#include "aslProbeTemplates.h"


namespace asl
{

	ProbeTemplate::ProbeTemplate (int n, AVec<double>* vec)
	{
		vectors.resize(n);
		for (int i(0); i < n; ++i)
			vectors[i] = vec[i];
	}
		
	AVec<double> probeD1q2Data[2] = {makeAVec(.5),makeAVec(-.5)};
	ProbeTemplate probeD1q2(2, probeD1q2Data);

	AVec<double> probeD2q3Data[3] = {makeAVec(0.,   1./sqrt(3.)),
					 makeAVec(-.5,-1./sqrt(3.)/2.),
					 makeAVec(.5, -1./sqrt(3.)/2.)};
	ProbeTemplate probeD2q3(3, probeD2q3Data);

	AVec<double> probeD3q4Data[4] = {makeAVec(0.,   0., sqrt(3.)/sqrt(8)),
					 makeAVec( 1./sqrt(3.),     0., -1./sqrt(24.)),
					 makeAVec(-1./sqrt(3.)/2.,-.5, -1./sqrt(24)),
					 makeAVec(-1./sqrt(3.)/2., .5, -1./sqrt(24))};
	ProbeTemplate probeD3q4(4, probeD3q4Data);


		
} // asl
