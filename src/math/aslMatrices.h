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


/// \file aslMatrices.h Matrices

#ifndef ASLMATRICES
#define ASLMATRICES


#include "../aslUtilities.h"
#include "aslVectors.h"


namespace asl
{
	/// class algebraic matrix. 
	/// The class is an implementation of a dynamic matrix with defined algebraic operations
	template <typename T = double> class AMatr
	{
		private: 
			unsigned int nRow;
			unsigned int nCol;			
			AVec<T> internalVec;
		public:
			inline AMatr();
			inline AMatr(unsigned int nR, unsigned int nC);
      		inline AMatr(const AMatr<T> &a);
			inline AMatr(unsigned int nR, unsigned int nC, AVec<T> v);
			template <typename T1> AMatr(const AMatr<T1> &a);
			const AMatr<T> & operator=(const AMatr & a);
			///doesn't chek boundaries
			inline T& operator()(int i, int j) {return internalVec[i*nCol+j];}
      		///doesn't chek boundaries
			inline const T& operator()(int i, int j)const {return internalVec[i*nCol+j];}
			///doesn't chek boundaries
			inline T& operator[](int i) {return internalVec[i];}
			///doesn't chek boundaries
			inline const T& operator[](int i)const {return internalVec[i];}
			inline unsigned int getNRow() const;
			inline unsigned int getNCol() const;
			inline void resize(unsigned int nR, unsigned int nCol);
			inline const AVec<T> & getInternalVec() const;
			inline AVec<T> & getInternalVec();
			void setRow(unsigned int r, const AVec<T> & a);
			void setColumn(unsigned int c, const AVec<T> & a);
	};

	
	/// \relates AMatr	
	template <typename T> std::ostream& operator<<(std::ostream &f,const AMatr<T> & a); 
	/// \relates AMatr	
	template <typename T> 
		inline const AMatr<T> & operator+=(AMatr<T> & a, const AMatr<T> & b);
	/// \relates AMatr	
	template <typename T> 
		inline const AMatr<T> operator+ (const AMatr<T> & a, const AMatr<T> & b); 
	/// \relates AMatr	
	template <typename T> 
		inline const AMatr<T> operator- (const AMatr<T> & a, const AMatr<T> & b);
	/// \relates AMatr	
	template <typename T> 
		const AMatr<T> operator* (const AMatr<T> &a, const AMatr<T> & b); 

	/// \relates AMatr	
	template <typename T> 
		const AVec<T> operator* (const AMatr<T> &a, const AVec<T> & b); 

	/// \relates AMatr	
	template <typename T> 
		const AVec<T> operator* (const AVec<T> &a, const AMatr<T> & b); 

	/// \relates AMatr	
	template <typename T> 
		const AMatr<T> operator* (const AMatr<T> &a, const T & b); 

	/// \relates AMatr	
	template <typename T> 
		const AMatr<T> operator* (const T &a, const AMatr<T> & b); 
	
	///Trace of a matrix \f$Tr(A)\equiv A_{ii}\f$ \relates AMatr	
	template <typename T> const T trace(const AMatr<T> &a); 

	///Trace of a matrix product \f$Tr(A B)\equiv A_{ij}B_{ji}\f$ \relates AMatr	
 	template <typename T> const T trace(const AMatr<T> & a, const AMatr<T> & b);
	/// \relates AMatr	
	template <typename T> 
		inline const AMatr<T> operator/ (const AMatr<T> & b, const T & a); 

	///element product of two vectors
	/**
		 \relates AMatr
		 \f$ elementProduct\left(
		                         \left[\begin{array}{c}
		                               a_1\\ \vdots \\ a_n
		                               \end{array}\right],
		                         \left[\begin{array}{c}
		                               b_1\\ \vdots \\ b_n
		                               \end{array}\right] =
		                         \left[\begin{array}{ccc}
		                               a_1b_1 & \cdots & a_1b_n\\
		                               \vdots & \ddots & \vdots\\
		                               a_nb_1 & \cdots & a_nb_n\\
		                               \end{array}\right]
		                         \right)
		 \f$
	 */			
	template <typename T>	
		AMatr<T> elementProduct(const AVec<T> & a, const AVec<T> & b);


	///	 generates a matrix with a row \relates AMatr
	template <typename T> AMatr<T> makeAMatr(const AVec<T> & a);
	
	///	 generates a matrix with two rows 	 \relates AMatr
	template <typename T> AMatr<T> makeAMatr(const AVec<T> & a, const AVec<T> & b);
	
	///	 generates a matrix with three rows  \relates AMatr
	template <typename T> AMatr<T> makeAMatr(const AVec<T> & a, 
	                                         const AVec<T> & b,
	                                         const AVec<T> & c);

	///	 generates a matrix with \p n rows  \relates AMatr
	template <typename T> AMatr<T> makeAMatr(AVec<T> *a, unsigned int n);

	/// \relates AMatr	
	template <typename T=int>AMatr<T> makeAMatrUnit(unsigned int n);

	
	/// returns AVec containing the diagonal elements
	/**
		\relates AMatr 
		the finction is valid only for square matrices
	*/
	template <typename T> AVec<T> getDiagonal(const AMatr<T> & a);
	
	/// returns AVec<T> containing the uper off diagonal elements
	/**
		\relates AMatr
		the function is valid only for square matrices
		\todo implement
	*/
	template <typename T> AVec<T> getOffDiagonalUp(const AMatr<T> & a);

	/// computes determinant expression fo cases 2x2 and 3x3 only
	/**
		 \relates AMatr		 
	*/
	template <typename T> T det(const AMatr<T> & m);

	/// returns solution of a system of linear equations
	/**
		\ingroup ComplexDataTypes		 
	*/
	template <typename T> AVec<T> solveSystem(const AMatr<T> & a, 
	                                          const AVec<T> & b);

	/// generate matrix with content of the matrix \p a but with replaced row \p r by vector \p b \relates AMatr<T>		 
	template <typename T> 
		AMatr<T> replaceRow(const AMatr<T> & a, const AVec<T> & b, unsigned int r);

	/// generate matrix with content of the matrix \p a but with replaced column \p c by vector \p b \relates AMatr<T>		 
	template <typename T> 
		AMatr<T> replaceColumn(const AMatr<T> & a, const AVec<T> & b, unsigned int c);

	/// returns inverse matrix for cases 2x2 and 3x3 \relates AMatr<T> 
	template <typename T> 
		AMatr<T> inverseMatrix(const AMatr<T> & a);

	
// ------------------------------ Implementation -----------------
	template <typename T> inline AMatr<T>::AMatr():
		nRow(0),
		nCol(0)
	{
	}

	template <typename T> inline AMatr<T>::AMatr(unsigned int nR, unsigned int nC):
		nRow(nR),
		nCol(nC),
		internalVec(nR*nC)
	{
	}
    template <typename T> inline AMatr<T>::AMatr(const AMatr<T> &a):
		nRow(a.nRow),
		nCol(a.nCol),
		internalVec(a.internalVec)
	{
	}
	template <typename T> inline AMatr<T>::AMatr(unsigned int nR, unsigned int nC, AVec<T> v):
		nRow(nR),
		nCol(nC),
		internalVec(v)
	{
	}

	template <typename T> inline AVec<T> & AMatr<T>::getInternalVec()
	{
		return internalVec;
	}

	template <typename T> inline const AVec<T> & AMatr<T>::getInternalVec() const
	{
		return internalVec;
	}
		
	template <typename T> inline unsigned int AMatr<T>::getNRow() const
	{
		return nRow;
	}
	
	template <typename T> inline unsigned int AMatr<T>::getNCol() const
	{
		return nCol;
	}
	
	template <typename T> 
		inline void AMatr<T>::resize(unsigned int nr, unsigned int nc)
	{
		nRow=nr; 
		nCol=nc;
		internalVec.resize(nr*nc);
	}

	template <typename T> 
		inline const AMatr<T> & operator+=(AMatr<T> & a, const AMatr<T> & b)
	{
		a.getInternalVec += b.getInternalVec;
		return a;
	}

	template <typename T> 
		inline const AMatr<T> operator+ (const AMatr<T> & a, const AMatr<T> & b) 
	{
		return {a.getNRow(),a.getNCol(),a.getInternalVec() + b.getInternalVec()};
	}

	template <typename T> 
		inline const AMatr<T> operator- (const AMatr<T> & a,const AMatr<T> & b) 
	{  
		return {a.getNRow(),a.getNCol(),a.getInternalVec() - b.getInternalVec()};
	}

	template <typename T> 
		inline const AMatr<T> operator/ (const AMatr<T> & b, const T & a) 
	{
		return {b.getNRow(), b.getNCol(), b.getInternalVec() / a};
	}
	
	/// Eigenvalues and eigenvectors calcutaion for symetric matrix 2x2
	/**
		\param a matrix element
		\param b matrix element
		\param c matrix element
		\param l1 first eigenvalue
		\param l2 second eigenvalue
		\param v1x x component of first eigenvector
		\param v1y y component of first eigenvector
		\param v2x x component of second eigenvector
		\param v2y y component of second eigenvector
		
		\f[
		     A= 
		     \left| \begin{array}{cc}
			a & c \\
			c & b 
		     \end{array}\right| 
		\f]
			
	*/
	void getEValEVecMatSym2x2(double a, double b, double c, double & l1, double & l2, double & v1x, double & v1y, double & v2x, double & v2y);
	/// Eigenvalues and eigenvectors calcutaion for symetric matrix 2x2
	/**
		\param a matrix element
		\param b matrix element
		\param c matrix element
		\param l1 first eigenvalue
		\param l2 second eigenvalue
		\param v1x x component of first eigenvector
		\param v1y y component of first eigenvector
		\param v2x x component of second eigenvector
		\param v2y y component of second eigenvector
		
		\f[
		     A= 
		     \left| \begin{array}{cc}
			a & c \\
			c & b 
		     \end{array}\right| 
		\f]
			
	*/
	void getEValEVecMatSym2x2(double a, double b, double c, double & l1, double & l2, double & v1x, double & v1y, double & v2x, double & v2y);

	/// Eigenvalues and eigenvectors calcutaion for symetric matrix 2x2
	/**
		\param a matrix element
		\param b matrix element
		\param c matrix element
		\param d matrix element
		\param e matrix element
		\param f matrix element
		\param l1 first eigenvalue
		\param l2 second eigenvalue
		\param l3 second eigenvalue
		\param v1x x component of first eigenvector
		\param v1y y component of first eigenvector
		\param v1z z component of first eigenvector
		\param v2x x component of second eigenvector
		\param v2y y component of second eigenvector
		\param v2z z component of second eigenvector
		\param v3x x component of second eigenvector
		\param v3y y component of second eigenvector
		\param v3z z component of second eigenvector
		
		\f[
		     A= 
		     \left| \begin{array}{ccc}
			a & d & e\\
			d & b & f \\
		    e & f & c \\
		     \end{array}\right| 
		\f]
			
	*/
	void getEValEVecMatSym3x3(double a, double b, double c, double e, double f, double g, 
	                          double & l1, double & l2, double & l3, 
	                          double & v1x, double & v1y, double & v1z, 
	                          double & v2x, double & v2y, double & v2z,
	                          double & v3x, double & v3y, double & v3z);
		 
} // asl

#endif //ASLMATRICES

