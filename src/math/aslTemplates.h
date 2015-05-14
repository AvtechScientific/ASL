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


#ifndef ASLTEMPLATES_H
#define ASLTEMPLATES_H

#include "aslVectors.h"
//#include "aslSVectors.h"

namespace asl {

/**
 \defgroup Templates Vector Templates
 */
	
/**
 \defgroup TemplatesNN Vector Templates: Nearest Neighbours
 \ingroup Templates
 */

/**
 \defgroup TemplatesNNP Vector Templates: Nearest Neighbours Plus
 \ingroup Templates
 */
	
/**
 \defgroup TemplatesEC Vector Templates: Elementary Cells
 \ingroup Templates
 */
	
/**
 \defgroup TemplatesNN0 Vector Templates: Nearest Neighbours without center
 \ingroup Templates
 */

/**
 \defgroup TemplatesNNP0 Vector Templates: Nearest Neighbours Plus without center
 \ingroup Templates
 */

/// list of implemented names of VectorTemplate
	enum VTName 
	{
		VTN_D1Q2EC,
		VTN_D2Q4EC,	
		VTN_D3Q8EC,
		VTN_D1Q3,	
		VTN_D2Q5,	
		VTN_D2Q9,	
		VTN_D3Q7,	
		VTN_D3Q15,	
		VTN_D3Q19,
		VTN_D3Q27,
		VTN_D1Q1UV,
		VTN_D2Q2UV,	
		VTN_D3Q3UV,
		VTN_D1Q2,
		VTN_D2Q4,
		VTN_D3Q6,
		VTN_D2Q8,
		VTN_D3Q14,
		VTN_D3Q18			
	};

	
	/// Defines set of vectros with several properties \ingroup Templates
	class VectorTemplate
	{
		private:
			void buildInvertVectorList();
		public:
			std::vector<AVec<int> > vectors;
			std::vector<double> laplasCoefs;
			std::vector<double> gradientCoefs;
			std::vector<double> quasiparticlesCoefs;
			std::vector<unsigned int> invertVectors;
			double dIdJLapCoef;
			std::vector<std::vector<double>> dxCoefs;
			std::vector<std::vector<std::vector<double>>> dIdJCoefs;
			VectorTemplate(int n, AVec<int>* vec);		
			VectorTemplate(int n, AVec<int>* vec, double* lc, double* gc);		

			inline unsigned int numberOfDimentions() const;
			inline AVec<int> getInverVector(unsigned int i);
	};

	inline unsigned int nD(const VectorTemplate & vt);
		
	///An elementary cell in 1D space
	/**	
	 \ingroup TemplatesEC
		 
	 This template defines an elementary cell: vectors and interpolation law
	 \image html t2ec.png "1D elementary cell"
	*/
	const VectorTemplate & d1q2ec();
	///An elementary cell in 2D space
	/**	
	 \ingroup TemplatesEC
		 
	 This template defines an elementary cell: vectors and interpolation law
	 \image html t4ec.png "2D elementary cell"
	*/
	const VectorTemplate & d2q4ec();	
	
	///An elementary cell in 3D space
	/**	
	 \ingroup TemplatesEC
		 
	 This template defines an elementary cell: vectors and interpolation law
	 \image html t8ec.png "3D elementary cell"
	*/
	const VectorTemplate & d3q8ec();
	
	///Vector template
	/**	
	 \ingroup TemplatesNN
		 
	 \image html t3.png
	*/	
	const VectorTemplate & d1q3();	

	///Vector template
	/**	
	 \ingroup TemplatesNN
		 
	 \image html t5.png
	*/	
	const VectorTemplate & d2q5();	

	///Vector template
	/**	
	 \ingroup TemplatesNNP
		 
	 \image html t9.png
	*/	
	const VectorTemplate & d2q9();	

	///Vector template
	/**	
	 \ingroup TemplatesNN
		 
	 \image html t7.png
	*/
	const VectorTemplate & d3q7();	

	///Vector template
	/**	
	 \ingroup TemplatesNNP
		 
	 \image html t15.png
	*/
	const VectorTemplate & d3q15();	

	///Vector template
	/**	
	 \ingroup TemplatesNNP
		 
	 \image html t19.png
	*/
	const VectorTemplate & d3q19();

	///Vector template
	/**	
	 \ingroup TemplatesNNP
		 
	 \image html t27.png
	*/
	const VectorTemplate & d3q27();
	
	///An unit vector in 1D space
	/**	
	 \ingroup TemplatesNNP
		 
	 This template defines a unit vector
	*/
	const VectorTemplate & d1q1uv();
	///An elementary cell in 2D space
	/**	
	 \ingroup Templates
		 
	 This template defines unit vectors
	*/
	const  VectorTemplate & d2q2uv();	
	
	///An elementary cell in 3D space
	/**	
	 \ingroup Templates
		 
	 This template defines unit vectors
	*/
	const VectorTemplate & d3q3uv();

	///Vector template
	/**	
	 \ingroup TemplatesNN0
		 
	 \image html t2.png
	*/	
	const VectorTemplate & d1q2();	

	///Vector template
	/**	
	 \ingroup TemplatesNN0
		 
	 \image html t4.png
	*/	
	const VectorTemplate & d2q4();	

	///Vector template
	/**	
	 \ingroup TemplatesNN0
		 
	 \image html t6.png
	*/	
	const VectorTemplate & d3q6();
	
	///Vector template
	/**	
	 \ingroup TemplatesNNP0
	*/	
	const VectorTemplate & d2q8();	

	///Vector template
	/**	
	 \ingroup TemplatesNNP0
 	*/	
	const VectorTemplate & d3q14();	
	
	///Vector template
	/**	
	 \ingroup TemplatesNNP0		 
	*/	
	const VectorTemplate & d3q18();	

	/// returns template corresponding to nearest neighbours \ingroup Templates
	/** asl::d1q3, asl::d2q5, asl::d3q7
	*/
	inline const VectorTemplate* nearestNeigboursVT(unsigned int dimNumber);

	/// returns template corresponding to nearest neighbours without central point
	/** 
		 \ingroup Templates
		 asl::d1q2, asl::d2q4, asl::d3q6
	*/
	inline const VectorTemplate* nearestNeigboursVT0(unsigned int dimNumber);

	/// returns template corresponding to nearest neighbours plus \ingroup Templates
	/** asl::d1q3, asl::d2q9, asl::d3q15
	*/
	inline const VectorTemplate* nearestNeigboursPVT(unsigned int dimNumber);

	
	/// returns template corresponding to all neighbours  
	/** 
		 \ingroup Templates
		 asl::d1q3, asl::d2q9, asl::d3q27
	*/
	inline const VectorTemplate* allNeigboursVT(unsigned int dimNumber);

	/// returns template corresponding to an elementary cell   
	/** 
		 \ingroup Templates
		 asl::d1q2ec, asl::d2q4ec, asl::d3q8ec
	*/
	inline const VectorTemplate* elementaryCellVT(unsigned int dimNumber);

	
// ----------------------------- Implementation -------------------------

	inline unsigned int VectorTemplate::numberOfDimentions() const
	{
		return vectors[0].getSize();
	}

	inline unsigned int nD(const VectorTemplate & vt)
	{
		return vt.numberOfDimentions();
	}

	
	inline const VectorTemplate* nearestNeigboursVT(unsigned int dimNumber)
	{
		static const VectorTemplate* vt[3]={&d1q3(),&d2q5(),&d3q7()};
		return vt[dimNumber-1];
	}

	inline const VectorTemplate* nearestNeigboursPVT(unsigned int dimNumber)
	{
		static const VectorTemplate* vt[3]={&d1q3(),&d2q9(),&d3q15()};
		return vt[dimNumber-1];
	}
	
	inline const VectorTemplate* allNeigboursVT(unsigned int dimNumber)
	{
		static const VectorTemplate* vt[3]={&d1q3(),&d2q9(),&d3q27()};
		return vt[dimNumber-1];
	}

	inline const VectorTemplate* elementaryCellVT(unsigned int dimNumber)
	{
		static const VectorTemplate* vt[3]={&d1q2ec(),&d2q4ec(),&d3q8ec()};
		return vt[dimNumber-1];
	}
	
	inline AVec<int> VectorTemplate::getInverVector(unsigned int i)
	{
		return vectors[invertVectors[i]];
	}
	
/*	
	///The templates for the numerical schemas
	namespace templ {
		/// An abstract class for description of point templates for differential operators etc.
		class GrigTemplate{
			public:
				unsigned const int nDimentions;
				unsigned const int nPoints;
				virtual vector<int> & getVector(unsigned int i)=0;
			protected:
				GridTemplate(unsigned int nD, unsigned int nP);	
		}
*/
		/**
		This template defines an elementary cell: vectors and interpolation law
		\image html t2ec.png "1D elementary cell"
		\image latex t2ec.png "1D elementary cell" width=8cm
		*/
/*		class d1q2ec {
			public:
				static const int ND=1;
				static const int NV=2;
				typedef Vec<ND,int> DV;
			public:
				inline static const DV & l(int i) {
					static const DV lat[NV] = {DV(0),DV(1)};
					return lat[i];
				}
				///The interpolation procedure
				template <typename T>
				inline static const T interpol(const Vec<NV,T> &d, const Vec<ND,lFl> &e) {
					return d.x()*(1.-e.x())+d.y()*e.x();
				}
		};
*/
		///An elementary cell in 2D space
		/**
		This template defines an elementary cell: vectors and interpolation law
		\image html t4ec.png "2D elementary cell"
		\image latex t4ec.png "2D elementary cell" width=8cm
		*/
/*		class d2q4ec {
			public:
				static const int ND=2;
				static const int NV=4;
				typedef Vec<ND,int> DV;
				inline static const DV & l(int i) {
					static const DV lat[NV] = {DV(0,0),DV(1,0),DV(0,1),DV(1,1)};
					return lat[i];
				}
				///The interpolation procedure
				template <typename T>
				inline static const T interpol(const Vec<NV,T> &d, const Vec<ND,lFl> &e) {
					return d[0]*(1.-e.x())*(1.-e.y())+d[1]*e.x()*(1.-e.y())+
					       d[2]*(1.-e.x())*e.y()     +d[3]*e.x()*e.y();
				}

		};
*/
		///An elementary cell in 2D space
		/**
		This template defines an elementary cell: vectors and interpolation law
		\image html t8ec.png "3D elementary cell"
		\image latex t8ec.png "3D elementary cell" width=8cm
		*/
/*		class d3q8ec {
			public:
				static const int ND=3;
				static const int NV=8;
				typedef Vec<ND,int> DV;
				inline static const DV & l(int i) {
					static const DV lat[NV] = {
							       DV(0,0,0),DV(1,0,0),DV(0,1,0),DV(1,1,0),
							       DV(0,0,1),DV(1,0,1),DV(0,1,1),DV(1,1,1)
							     };
					return lat[i];
				}

				///The interpolation procedure
				template <typename T>
				inline static const T interpol(const Vec<NV,T> &d, const Vec<ND,lFl> &e) {
					return d[0]*(1.-e.x())*(1.-e.y())*(1.-e.z())
					 +d[1]*    e.x()* (1.-e.y())*(1.-e.z())
					 +d[2]*(1.-e.x())*    e.y()* (1.-e.z())
					 +d[3]*    e.x()*     e.y()* (1.-e.z())
					 +d[4]*(1.-e.x())*(1.-e.y())*    e.z()
					 +d[5]*    e.x()* (1.-e.y())*    e.z()
					 +d[6]*(1.-e.x())*    e.y()*     e.z()
					 +d[7]*    e.x()*     e.y()*     e.z();
				}
		};

*/

		///d2q9 vectorspace
		/**
		The operators are defined on the correspoNDing templates:
		\image html t9.png "9-point templates"
		\image latex t9.png "9-point templates" width=8cm
		*/
/*		class d2q9 {
			public:
				static const int ND=2;
				static const int NV=9;
				typedef Vec<ND,int> DV;
				//      private:
				static const lFl w1=4./9.,w2=1./9.,w3=1./36., as2i=3.;
				static inline const DV & l(int i) {
					static const DV lat[NV] = {
							       DV(0,0),DV(1,0),DV(0,1), DV(-1,0),DV(0,-1),
							       DV(1,1),DV(-1,1),DV(-1,-1),DV(1,-1)
							     };
					return lat[i];
				}
			///number of the vector wich has oposite derection
				static inline int io(int i) {
				  static const int inv[NV] = {0,3,4,1,2,7,8,5,6};
				  return inv[i];
				}
				static inline lFl s(const Vec<NV> &p) {
					return (w1*p[0]+w2*(p[1]+p[2]+p[3]+p[4])+w3*(p[5]+p[6]+p[7]+p[8]));
				}
				static inline lFl s(Vec<NV> &p, const lFl &val) {
					return (p[0] =p[1] =p[2] =p[3] =p[4] =p[5] =p[6] =p[7] =p[8] =val);
				}

				static inline const Vec<ND> v(const Vec<NV> &p) {
				  return Vec<2>(w2*(p[1]-p[3])+w3*(p[5]-p[6]-p[7]+p[8]),
						w2*(p[2]-p[4])+w3*(p[5]+p[6]-p[7]-p[8]));
				}
				static inline const Vec<ND> v(Vec<NV> &p, const Vec<ND> &val) {
					const lFl p0(s(p));
					for (int i(0); i < NV; ++i)
						p[i] =p0+val*Vec<ND>(l(i))*as2i;
					return val;
				}

				///coefficient of the lapplace operator
				static inline lFl lc(int i) {
					static const lFl c[NV] = {-10./3.,2./3.,2./3.,2./3.,2./3.,1./6.,1./6.,1./6.,1./6.};
					return c[i];
				}

		};
*/

		///d2q5 vectorspace
		/**
		The operators are defined on the correspoNDing templates:
		\image html t5.png "5-point templates"
		\image latex t5.png "5-point templates" width=8cm
		*/
/*		class d2q5 {
			public:
				static const int ND=2;
				static const int NV=5;
				typedef Vec<ND,int> DV;
				inline static const DV & l(int i) {
					static const DV lat[NV] = {DV(0,0),DV(1,0),DV(0,1),DV(-1,0), DV(0,-1)};
					return lat[i];
				}
				///number of the vector wich has oposite derection
				inline static int io(int i) {
					static const int inv[NV] = {0,3,4,1,2};
					return inv[i];
				}
				///coefficient of the lapplace operator
				static inline lFl lc(int i) {
					static const lFl c[NV] = {-4.,1.,1.,1.,1.};
					return c[i];
				}
		};
*/

		///d2q4 vectorspace
		/**
		The operators are defined on the corresponding templates:
		\image html t4.png "4-point templates"
		\image latex t4.png "4-point templates" width=8cm
		*/
/*		class d2q4 {
			public:
				static const int ND=2;
				static const int NV=4;
				typedef Vec<ND,int> DV;
				static const lFl w1=1./NV, as2i=ND;
				inline static const DV & l(int i) {
					static const DV lat[NV] = {DV(1,0),DV(0,1),DV(-1,0), DV(0,-1)};
					return lat[i];
				}
				///number of the vector wich has oposite derection
				inline static int io(int i) {
					static const int inv[NV] = {2,3,0,1};
					return inv[i];
				}
				static inline lFl s(const Vec<NV> &p) {return (w1*(p[0]+p[1]+p[2]+p[3]));}
				static inline lFl s(Vec<NV> &p, const lFl &val) {return (p[0] =p[1] =p[2] =p[3] =val);}
				static inline const Vec<ND> v(const Vec<NV> &p) {
					return Vec<ND>(w1*(p[0]-p[2]),w1*(p[1]-p[3]));
				}
				static inline const Vec<ND> v(Vec<NV> &p, const Vec<ND> &val) {
					const lFl p0(s(p));
					for (int i(0); i < NV; ++i)p[i] =p0+val*Vec<ND>(l(i))*as2i;
					return val;
				}

		};

*/

		///d3q7 vectorspace
		/**
		\image html t7.png "7-point templates"
		\image latex t7.png "7-point templates" width=8cm
		*/
/*		class d3q7 {
			public:
				static const int ND=3;
				static const int NV=7;
				typedef Vec<ND,int> DV;
				inline static const DV l(int i) {
					static const DV lat[NV] = {DV(SV()[I2T<0>()]),DV(SV()[I2T<1>()]),DV(SV()[I2T<2>()]),
						      DV(SV()[I2T<3>()]),DV(SV()[I2T<4>()]),DV(SV()[I2T<5>()]),
						      DV(SV()[I2T<6>()])};
					return lat[i];
				}
				///number of the vector wich has oposite derection
				inline static int io(int i) {
					static const int inv[NV] = {0,4,5,6,1,2,3};
					return inv[i];
				}
				///coefficient of the lapplace operator
				static inline lFl lc(int i) {
					static const lFl c[NV] = {-6.,1.,1.,1.,1.,1.,1.};
					return c[i];
				}
		};
*/

		///d3q6 vectorspace
		/**
		\image html t6.png "6-point templates"
		\image latex t6.png "6-point templates" width=8cm
		*/
/*		class d3q6 {
			public:
				static const int ND=3;
				static const int NV=6;
				typedef Vec<ND,int> DV;
				static const lFl w1=1./NV, as2i=ND;
				inline static const DV & l(int i) {
					static const DV lat[NV] = {
							       DV(1,0,0),   DV(0,1,0),DV(0,0,1),
							       DV(-1,0,0),  DV(0,-1,0),DV(0,0,-1)
							     };
					return lat[i];
				}
				///number of the vector wich has oposite derection
				inline static int io(int i) {
					static const int inv[NV] = {3,4,5,0,1,2};
					return inv[i];
				}
				static inline lFl s(const Vec<NV> &p) {return (w1*(p[0]+p[1]+p[2]+p[3]+p[4]+p[5]));}
				static inline lFl s(Vec<NV> &p, const lFl &val) {
					return (p[0] =p[1] =p[2] =p[3] =p[4] =p[5] =val);
				}
				static inline const Vec<ND> v(const Vec<NV> &p) {
					return Vec<ND>(w1*(p[0]-p[3]),w1*(p[1]-p[4]),w1*(p[2]-p[5]));
				}
				static inline const Vec<ND> v(Vec<NV> &p, const Vec<ND> &val) {
					const lFl p0(s(p));
					for (int i(0); i < NV; ++i)p[i] =p0+val*Vec<ND>(l(i))*as2i;
					return val;
				}
		};

*/

		///d3q19 vectorspace
		/**
		\image html t19.png "19-point templates"
		\image latex t19.png "19-point templates" width=8cm
		*/
/*		class d3q19 {
			public:
				static const int ND=3;
				static const int NV=19;
				typedef Vec<ND,int> DV;
				static const lFl w1=1./3.,w2=1./18.,w3=1./36., as2i=3.;
				inline static const DV & l(int i) {
					static const DV lat[NV] = {
							       DV(0,0,0),   DV(1,0,0),   DV(0,1,0),
							       DV(0,0,1),   DV(-1,0,0),  DV(0,-1,0),
							       DV(0,0,-1),  DV(1,1,0),   DV(1,-1,0),
							       DV(-1,-1,0), DV(-1,1,0),  DV(0,1,1),
							       DV(0,1,-1),  DV(0,-1,-1), DV(0,-1,1),
							       DV(1,0,1),   DV(1,0,-1),  DV(-1,0,-1),
							       DV(-1,0,1)
							     };
					return lat[i];
				}
				///numbers of the vector wich has oposite derection
				inline static int io(int i) {
					static const int inv[NV] = {0,4,5,6,1,2,3,9,10,7,8,13,14,11,12,17,18,15,16};
					return inv[i];
				}
				static inline lFl s(const Vec<NV> &p) {
					return (w1*p[0]+w2*(p[1]+p[2]+p[3]+p[4]+p[5]+p[6])
					+w3*(p[7]+p[8]+p[9]+p[10]+p[11]+p[12]+p[13]+p[14]+p[15]+p[16]+p[17]+p[18]));
				}
				static inline lFl s(Vec<NV> &p, const lFl &val) {p=val; return val;}

				static inline const Vec<ND> v(const Vec<NV> &p) {
					return 
					Vec<ND>(w2*(p[1]-p[4])+w3*( p[7] +p[8]- p[9]-p[10]+p[15]+p[16]-p[17]-p[18]),
						w2*(p[2]-p[5])+w3*( p[7] -p[8]- p[9]+p[10]+p[11]+p[12]-p[13]-p[14]),
						w2*(p[3]-p[6])+w3*(p[11]-p[12]-p[13]+p[14]+p[15]-p[16]-p[17]+p[18]));
				}
				static inline const Vec<ND> v(Vec<NV> &p, const Vec<ND> &val) {
					const lFl p0(s(p));
					for (int i(0); i < NV; ++i)p[i] =p0+val*Vec<ND>(l(i))*as2i;
					return val;
				}
				///coefficient of the lapplace operator
				static inline lFl lc(int i) {
				  static const lFl c[NV] = {-4.,
							   1./3.,1./3.,1./3.,1./3.,1./3.,1./3.,
							   1./6.,1./6.,1./6.,1./6.,1./6.,1./6.,
							   1./6.,1./6.,1./6.,1./6.,1./6.,1./6.};
				  return c[i];
				}
		};
*/

		///d3q15 vectorspace
		/**
		\image html t15.png "15-point templates"
		\image latex t15.png "15-point templates" width=8cm
		*/
/*		class d3q15 {
			public:
				static const int ND=3;
				static const int NV=15;
				typedef Vec<ND,int> DV;
				static const lFl w1=2./9.,w2=1./9.,w3=1./72., as2i=3.;
				static inline const DV & l(int i) {
					static const DV lat[NV] = {
						       DV(0,0,0),   DV(1,0,0),   DV(0,1,0),
						       DV(0,0,1),   DV(-1,0,0),  DV(0,-1,0),
						       DV(0,0,-1),  DV(1,1,1),   DV(1,1,-1),
						       DV(1,-1,1),  DV(1,-1,-1), DV(-1,1,1),
						       DV(-1,1,-1), DV(-1,-1,1), DV(-1,-1,-1)
						     };
					return lat[i];
				}
				///number of the vector wich has oposite derection
				static inline int io(int i) {
					static const int inv[NV] = {0,4,5,6,1,2,3,14,13,12,11,10,9,8,7};
					return inv[i];
				}
				static inline lFl s(const Vec<NV> &p) {
					return (w1*p[0]+w2*(p[1]+p[2]+p[3]+p[4]+p[5]+p[6])
						+w3*(p[7]+p[8]+p[9]+p[10]+p[11]+p[12]+p[13]+p[14]));
				}
				static inline lFl s(Vec<NV> &p, const lFl &val) {p=val; return val;}

				template <typename Td,typename SD> 
				static inline const Vec<ND> v(const Vec<NV,Td,SD> &p) 
				{
					return 
					Vec<ND>(w2*(p[1]-p[4])+w3*(p[7]+p[8]+p[9]+p[10]-p[11]-p[12]-p[13]-p[14]),
						w2*(p[2]-p[5])+w3*(p[7]+p[8]-p[9]-p[10]+p[11]+p[12]-p[13]-p[14]),
						w2*(p[3]-p[6])+w3*(p[7]-p[8]+p[9]-p[10]+p[11]-p[12]+p[13]-p[14]));
				}
				static inline void v(Vec<NV> &p, const Vec<ND> &val) {
					const lFl p0(s(p));
					for (int i(0); i < NV; ++i)p[i] =p0+val*Vec<ND>(l(i))*as2i;
				}
				template <typename Ty> 
				static inline void v(Vec<NV,lFl,Ty> p, const Vec<ND> &val) 
				{
					const lFl p0(s(p));
					for (int i(0); i < NV; ++i)p[i] =p0+val*Vec<ND>(l(i))*as2i;
				}

				///coefficient of the lapplace operator
				static inline lFl lc(int i) {
					static const lFl c[NV] = {-14./3.,
								2./3.,2./3.,2./3.,2./3.,2./3.,2./3.,
								1./12.,1./12.,1./12.,1./12.,
								1./12.,1./12.,1./12.,1./12.};
					return c[i];
				}
		};
*/

		///D3Q27 vectorspace
		/**
		\image html t27.png "27-point templates"
		\image latex t27.png "27-point templates" width=8cm
		*/
/*		class d3q27 {
			public:
				static const int ND=3;
				static const int NV=27;
				typedef Vec<ND,int> DV;
				inline static const DV & l(int i) {
					static const DV lat[NV] = {
							       DV(0,0,0),   DV(1,0,0),   DV(0,1,0),
							       DV(0,0,1),   DV(-1,0,0),  DV(0,-1,0),
							       DV(0,0,-1),  DV(1,1,0),   DV(1,-1,0),
							       DV(-1,-1,0), DV(-1,1,0),  DV(0,1,1),
							       DV(0,1,-1),  DV(0,-1,-1), DV(0,-1,1),
							       DV(1,0,1),   DV(1,0,-1),  DV(-1,0,-1),
							       DV(-1,0,1),  DV(1,1,1),   DV(1,1,-1),
							       DV(1,-1,1),  DV(1,-1,-1), DV(-1,1,1),
							       DV(-1,1,-1), DV(-1,-1,1), DV(-1,-1,-1)
							     };
					return lat[i];
				}
				///number of the vector wich has oposite derection
				inline static int io(int i) {
					static const int inv[NV] = {0,
								4,5,6,1,2,3,
								9,10,7,8,13,14,11,12,17,18,15,16,
								26,25,24,23,22,21,20,19};
					return inv[i];
				}
		};
*/

		///D1Q3 vectorspace
		/**
		\image html t3.png "3-point templates"
		\image latex t3.png "3-point templates" width=8cm
		*/
/*		class d1q3 {
			public:
				static const int ND=1; ///<number of dimensions
				static const int NV=3; ///<nuber of directions
				typedef Vec<ND,int> DV;      ///< type of lattice vector
				inline static const DV & l(int i) { /// Returns value of the lattice vector
					static const DV lat[NV] = {DV(0),DV(1),DV(-1)};
					return lat[i];
				}
				///number of the vector wich has oposite derection
				inline static int io(int i) {
					static const int inv[NV] = {0,2,1};
					return inv[i];
				}
		};
*/

		///D1Q2 vectorspace
		/**
		\image html t2.png "2-point templates"
		\image latex t2.png "2-point templates" width=8cm
		*/
/*		class d1q2 {
			public:
				static const int ND=1; ///<number of dimensions
				static const int NV=2; ///<nuber of directions
				typedef Vec<1,int> DV;      ///< type of lattice vector
				inline static const DV & l(int i) { /// Returns value of the lattice vector
					static const DV lat[NV] = {DV(1),DV(-1)};
					return lat[i];
				}
				///number of the vector wich has oposite derection
				inline static int io(int i) {
					static const int inv[NV] = {1,0};
					return inv[i];
				}
		};

		///Reflected template vector \f$ lrf(f)_i=f_{\tilde i}, \f$ \f$ \vec a_{\tilde i}=-\vec a_i \f$
		template <typename Tl> inline Vec<Tl::NV> lrv(const Vec<Tl::NV> &f_) {
			Vec<Tl::NV> a;
			for (int i(0); i < Tl::NV; ++i) {a[i] =f_[Tl::io(i)];}
			return a;
		}

		template <typename Tl,int I> inline void _l1f_c(I2T<I>,Tl, const Vec<Tl::ND> &v_,Vec<Tl::NV> &a){
			typedef typename Tl::SV SV;
			a[I] =SV()[I2T<I>()]*v_; _l1f_c(I2T<I-1>(),Tl(),v_,a);
		}
		template <typename Tl> inline void _l1f_c(I2T<0>,Tl, const Vec<Tl::ND> &v_,Vec<Tl::NV> &a){
			typedef typename Tl::SV SV; 
			a[0] =SV()[I2T<0>()]*v_;
		}

		///Creates template vector as scalar product with lattice vectors \f$ \vec a_i \vec v \f$
		template <typename Tl> inline Vec<Tl::NV> l1f(const Vec<Tl::ND> &v_) {
			typedef Vec<Tl::ND> TV;
			Vec<Tl::NV> a;
			//      for (int i(0); i < Tl::NV; ++i) {a[i] =TV(Tl::l(i))*v_;}
			_l1f_c(I2T<Tl::NV-1>(),Tl(),v_,a);
			return a;
		}

		template <typename Tl,int I> 
		inline void _l2f_c(I2T<I>,Tl, const Vec<Tl::ND> &v1_, const Vec<Tl::ND> &v2_,Vec<Tl::NV> &a){
			class Tl::SV q; a[I] =(q[I2T<I>()]*v1_)*(q[I2T<I>()]*v2_); _l2f_c(I2T<I-1>(),Tl(),v1_,v2_,a);
		}
		template <typename Tl> 
		inline void _l2f_c(I2T<0>,Tl, const Vec<Tl::ND> &v1_, const Vec<Tl::ND> &v2_,Vec<Tl::NV> &a){
			class Tl::SV q; a[0] =(q[I2T<0>()]*v1_)*(q[I2T<0>()]*v2_);
		}

		///Creates template vector as 2 scalar products with lattice vectors \f$ (\vec a_i \vec v_1) (\vec a_i \vec v_2) \f$
		template <typename Tl> inline Vec<Tl::NV> l2f(const Vec<Tl::ND> &v1_, const Vec<Tl::ND> &v2_) {
			typedef Vec<Tl::ND> TV;
			Vec<Tl::NV> a;
			//      for (int i(0); i < Tl::NV; ++i) {a[i] =(TV(Tl::l(i))*v1_)*(TV(Tl::l(i))*v2_);}
			_l2f_c(I2T<Tl::NV-1>(),Tl(),v1_,v2_,a);
			return a;
		}


	} // templ
*/
}// asl

#endif // TEMPL_H_INCLUDED
