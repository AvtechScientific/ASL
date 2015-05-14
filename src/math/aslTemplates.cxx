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


#include "aslTemplates.h"


namespace asl
{
	
	VectorTemplate::VectorTemplate (int n, AVec<int>* vec)
	{
		vectors.resize(n);
		for (int i(0); i < n; ++i)
			vectors[i] = vec[i];
		buildInvertVectorList();
	}

	VectorTemplate::VectorTemplate (int n, AVec<int>* vec, double* lc, double* gc):
		dIdJLapCoef(0.5)
	{
		vectors.resize(n);
		laplasCoefs.resize(n);
		gradientCoefs.resize(n);
		quasiparticlesCoefs.resize(n);
		unsigned int nD(vec[0].getSize());
		dxCoefs.resize(nD);
		for (unsigned int i(0); i < nD; ++i)
			dxCoefs[i].resize(n);
		dIdJCoefs.resize(nD);
		for (unsigned int i(0); i < nD; ++i){
			dIdJCoefs[i].resize(nD);
			for (unsigned int j(0); j<nD; ++j)
				dIdJCoefs[i][j].resize(n);
		}
		for (int i(0); i < n; ++i)
		{
			vectors[i] = vec[i];
			laplasCoefs[i] = lc[i];
			gradientCoefs[i] = gc[i];
			quasiparticlesCoefs[i] = lc[i]/6.;
			for (unsigned int k(0); k < nD; ++k)
			{
				dxCoefs[k][i] = gradientCoefs[i] * vectors[i][k];
				for (unsigned int m(0); m < nD; ++m)
					dIdJCoefs[k][m][i] = gradientCoefs[i] * vectors[i][k] * vectors[i][m] * 3.;
			}
		}
		quasiparticlesCoefs[0] += 1.;
		buildInvertVectorList();
	}

	void VectorTemplate::buildInvertVectorList()
	{
		unsigned int n(vectors.size());
		invertVectors.resize(n);
		for(unsigned int i(0); i < n; ++i)
			for(unsigned int j(i+1); j < n; ++j)
				if(vectors[i] == -vectors[j])
				{
					invertVectors[i] = j;
					invertVectors[j] = i;
				}
	}
		
	const VectorTemplate & d1q2ec()
	{	
		static AVec<int> d1q2ecData[2] = {makeAVec(0),makeAVec(1)};
		static VectorTemplate vt(2, d1q2ecData);
		return  vt;
	}

	const VectorTemplate & d2q4ec()
	{
		static AVec<int> d2q4ecData[4] = {makeAVec(0,0),makeAVec(1,0),makeAVec(0,1),makeAVec(1,1)};
		static VectorTemplate vt(4, d2q4ecData);
		return vt;
	}	

	const VectorTemplate & d3q8ec()
	{
		static AVec<int> d3q8ecData[8] ={makeAVec(0,0,0),makeAVec(1,0,0),makeAVec(0,1,0),makeAVec(1,1,0),
							 makeAVec(0,0,1),makeAVec(1,0,1),makeAVec(0,1,1),makeAVec(1,1,1)};
		static VectorTemplate vt(8, d3q8ecData);
		return vt;
	}	

	const VectorTemplate & d1q3()
	{
		static AVec<int> d1q3DataV[3] ={makeAVec(0),makeAVec(1),makeAVec(-1)};
		static double d1q3DataLC[3] = {-2,1,1};	
		static double d1q3DataGC[3] = {0,.5,.5};	
		static VectorTemplate vt(3, d1q3DataV,d1q3DataLC,d1q3DataGC);
		return vt;
	}	

		
	const VectorTemplate & d2q5()
	{
		static AVec<int> d2q5DataV[5] ={makeAVec(0,0),makeAVec(1,0),makeAVec(0,1),makeAVec(-1,0), makeAVec(0,-1)};
		static double d2q5DataLC[5] = {-4,1,1,1,1};	
		static double d2q5DataGC[5] = {0,.5,.5,.5,.5};	
		static VectorTemplate vt(5, d2q5DataV,d2q5DataLC,d2q5DataGC);
		return vt;
	}	

	const VectorTemplate & d2q9()
	{
		static AVec<int> d2q9DataV[9] = {makeAVec(0,0),makeAVec(1,0),makeAVec( 0,1),makeAVec(-1, 0), makeAVec(0,-1),
			                                           makeAVec(1,1),makeAVec(-1,1),makeAVec(-1,-1), makeAVec(1,-1)};
		static double d2q9DataLC[9] = {-10./3.,2./3.,2./3.,2./3.,2./3.,1./6.,1./6.,1./6.,1./6.};	
		static double d2q9DataGC[9] = {0,1./3.,1./3.,1./3.,1./3.,1./12.,1./12.,1./12.,1./12.};	
		static VectorTemplate vt(9, d2q9DataV,d2q9DataLC,d2q9DataGC);
		return vt;
	}	
		
	const VectorTemplate & d3q7()
	{
		static AVec<int> d3q7DataV[7] = {makeAVec(0,0,0),makeAVec(1,0,0), makeAVec(0,1,0), makeAVec(0,0,1),
												  makeAVec(-1,0,0),makeAVec(0,-1,0),makeAVec(0,0,-1)};
		static double d3q7DataLC[7] = {-6,1,1,1,1,1,1};	
		static double d3q7DataGC[7] = {0,.5,.5,.5,.5,.5,.5};	
		static VectorTemplate vt(7, d3q7DataV,d3q7DataLC,d3q7DataGC);
		return vt;
	}	

	const VectorTemplate & d3q15()
	{
		static AVec<int> d3q15DataV[15] = {makeAVec(0,0,0),makeAVec(1,0,0),  makeAVec(0,1,0), makeAVec(0,0,1),
											  makeAVec(-1,0,0), makeAVec(0,-1,0),makeAVec(0,0,-1),
											  makeAVec(1,1,1),  makeAVec(1,1,-1),makeAVec(1,-1,1),
											  makeAVec(1,-1,-1),makeAVec(-1,1,1),makeAVec(-1,1,-1),
											  makeAVec(-1,-1,1),makeAVec(-1,-1,-1)};
		static double d3q15DataLC[15] = {-14./3.,2./3.,2./3.,2./3.,
										2./3.,2./3.,2./3.,
										1./12.,1./12.,1./12.,
										1./12.,1./12.,1./12.,
										1./12.,1./12.};	
		static double d3q15DataGC[15] = {0,	1./3.,1./3.,1./3.,
										1./3.,1./3.,1./3.,
										1./24.,1./24.,1./24.,
										1./24.,1./24.,1./24.,
										1./24.,1./24.};	
		static VectorTemplate vt(15, d3q15DataV,d3q15DataLC,d3q15DataGC);
		return vt;
	}	
	

	const VectorTemplate & d3q19()
	{
		static AVec<int> d3q19DataV[19] = {makeAVec(0,0,0),makeAVec(1,0,0),  makeAVec(0,1,0),  makeAVec(0,0,1),
											  makeAVec(-1,0,0), makeAVec(0,-1,0), makeAVec(0,0,-1),
											  makeAVec(1,1,0),  makeAVec(1,-1,0), makeAVec(-1,-1,0),
											  makeAVec(-1,1,0), makeAVec(0,1,1),  makeAVec(0,1,-1),
											  makeAVec(0,-1,-1),makeAVec(0,-1,1), makeAVec(1,0,1),
											  makeAVec(1,0,-1), makeAVec(-1,0,-1),makeAVec(-1,0,1)};

		static double d3q19DataLC[19] = {-4.,1./3.,1./3.,1./3.,1./3.,1./3.,1./3.,
									1./6.,1./6.,1./6.,1./6.,1./6.,1./6.,
									1./6.,1./6.,1./6.,1./6.,1./6.,1./6.};	
		static double d3q19DataGC[19] = {0,	1./6.,1./6.,1./6.,1./6.,1./6.,1./6.,
									1./12.,1./12.,1./12.,1./12.,1./12.,1./12.,
									1./12.,1./12.,1./12.,1./12.,1./12.,1./12.};	
		static VectorTemplate vt(19, d3q19DataV,d3q19DataLC,d3q19DataGC);
		return vt;
	}	

	const VectorTemplate & d3q27()
	{
		static AVec<int> d3q27DataV[27] = {makeAVec(0,0,0),makeAVec(1,0,0),  makeAVec(0,1,0),  makeAVec(0,0,1),
												  makeAVec(-1,0,0), makeAVec(0,-1,0), makeAVec(0,0,-1),
												  makeAVec(1,1,0),  makeAVec(1,-1,0), makeAVec(-1,-1,0),
												  makeAVec(-1,1,0), makeAVec(0,1,1),  makeAVec(0,1,-1),
												  makeAVec(0,-1,-1),makeAVec(0,-1,1), makeAVec(1,0,1),
												  makeAVec(1,0,-1), makeAVec(-1,0,-1),makeAVec(-1,0,1),
												  makeAVec(1,1,1),  makeAVec(1,1,-1),makeAVec(1,-1,1),
												  makeAVec(1,-1,-1),makeAVec(-1,1,1),makeAVec(-1,1,-1),
												  makeAVec(-1,-1,1),makeAVec(-1,-1,-1)};
		static VectorTemplate vt(27, d3q27DataV);
		return vt;
	}	
		
	const VectorTemplate & d1q1uv()
	{
		static AVec<int> d1q1uvData[1] = {makeAVec(1)};
		static VectorTemplate vt(1, d1q1uvData);
		return vt;
	}	

	const VectorTemplate & d2q2uv()
	{
		static AVec<int> d2q2uvData[2] = {makeAVec(1,0),makeAVec(0,1)};
		static VectorTemplate vt(2, d2q2uvData);
		return vt;
	}	

	const VectorTemplate & d3q3uv()
	{
		static AVec<int> d3q3uvData[3] ={makeAVec(1,0,0),makeAVec(0,1,0),makeAVec(0,0,1)};
		static VectorTemplate vt(3, d3q3uvData);
		return vt;
	}	

	const VectorTemplate & d1q2()
	{
		static AVec<int> d1q2DataV[2] ={makeAVec(1),makeAVec(-1)};	
		static VectorTemplate vt(2, d1q2DataV);
		return vt;
	}	

	const VectorTemplate & d2q4()
	{
		static AVec<int> d2q4DataV[4] ={makeAVec(1,0),makeAVec(0,1),makeAVec(-1,0), makeAVec(0,-1)};
		static VectorTemplate vt(4, d2q4DataV);
		return vt;
	}	
		
	const VectorTemplate & d3q6()
	{
		static AVec<int> d3q6DataV[6] = {makeAVec(1,0,0), makeAVec(0,1,0), makeAVec(0,0,1),
								makeAVec(-1,0,0),makeAVec(0,-1,0),makeAVec(0,0,-1)};
		static VectorTemplate vt(6, d3q6DataV);
		return vt;
	}	


	const VectorTemplate & d2q8()
	{
		static AVec<int> d2q8DataV[8] ={makeAVec(1,0),makeAVec(0,1),makeAVec(-1,0), makeAVec(0,-1),
							 makeAVec(1,1),makeAVec(-1,1),makeAVec(-1,-1),makeAVec(1,-1)};
		static VectorTemplate vt(8, d2q8DataV);
		return vt;
	}	
		
	const VectorTemplate & d3q14()
	{
		static AVec<int> d3q14DataV[14] = {makeAVec(1,0,0),  makeAVec(0,1,0), makeAVec(0,0,1),
								makeAVec(-1,0,0), makeAVec(0,-1,0),makeAVec(0,0,-1),
								makeAVec(1,1,1),  makeAVec(1,1,-1),makeAVec(1,-1,1),
								makeAVec(1,-1,-1),makeAVec(-1,1,1),makeAVec(-1,1,-1),
								makeAVec(-1,-1,1),makeAVec(-1,-1,-1)};
		static VectorTemplate vt(14, d3q14DataV);
		return vt;
	}	

	const VectorTemplate & d3q18()
	{
		static AVec<int> d3q18DataV[18] = {makeAVec(1,0,0),  makeAVec(0,1,0),  makeAVec(0,0,1),
							    makeAVec(-1,0,0), makeAVec(0,-1,0), makeAVec(0,0,-1),
								makeAVec(1,1,0),  makeAVec(1,-1,0), makeAVec(-1,-1,0),
								makeAVec(-1,1,0), makeAVec(0,1,1),  makeAVec(0,1,-1),
								makeAVec(0,-1,-1),makeAVec(0,-1,1), makeAVec(1,0,1),
								makeAVec(1,0,-1), makeAVec(-1,0,-1),makeAVec(-1,0,1)};
		static VectorTemplate vt(18, d3q18DataV);
		return vt;
	}	
		
}// asl
