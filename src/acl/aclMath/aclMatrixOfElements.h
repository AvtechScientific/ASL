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


#ifndef ACLMATRIXOFELEMENTS_H
#define ACLMATRIXOFELEMENTS_H

#include "aclVectorOfElementsDef.h"

namespace acl
{
	/// The class represents a matrix elements of ::Element
	/**
		 \ingroup ComplexDataTypes
	 */
	class MatrixOfElements
	{
		private:
			/// number of columns 
			unsigned int nRow;
			/// number of rows
			unsigned int nCol;

			VectorOfElements ve; 
			/// defines convertion rule of matrix indeces \p i and \p j into vector one
			unsigned int ij2i(unsigned int i, unsigned int j) const;
		public:	
			explicit MatrixOfElements(unsigned int nR = 0, unsigned int nC = 0);
			
			void setElement(unsigned int r, unsigned int c, Element a);
			void setRow(unsigned int r,const VectorOfElements & a);
			void setColumn(unsigned int c, const VectorOfElements & a);			
			const Element getElement(unsigned int r, unsigned int c) const;
			const VectorOfElements getVE(unsigned int r, unsigned int c) const;
			const unsigned int getNColumns() const;
			const unsigned int getNRows() const;
			VectorOfElements & getInternalVector();
			const VectorOfElements & getInternalVector() const;
			inline void resize(unsigned int nr, unsigned int nc);
			MatrixOfElements operator= (const MatrixOfElements & m);
	};

	///function copies the MatrixOfElements class. 
	/**
		 \relates MatrixOfElements
	 */
	void copy(const MatrixOfElements & source, MatrixOfElements & destination);

	///summ of two matrices
	/**
		 \relates MatrixOfElements
	 */	
	MatrixOfElements operator+(const MatrixOfElements & a, const MatrixOfElements & b);
	///difference of two matrices
	/**
		 \relates MatrixOfElements
	 */	
	MatrixOfElements operator-(const MatrixOfElements & a, const MatrixOfElements & b);	

	///product of two matrices
	/**
		 \relates MatrixOfElements
	 */	
	MatrixOfElements operator*(const MatrixOfElements & a, const MatrixOfElements & b);	

	///product of vector and matrix
	/**
		 \relates MatrixOfElements
	 */		
	VectorOfElements operator*(const VectorOfElements & a, const MatrixOfElements & b);

	///product of vector and matrix
	/**
		 \relates MatrixOfElements
	 */		
	VectorOfElements operator*(const MatrixOfElements & a, const VectorOfElements & b);

	///division of a matrix on a VectorOfElements with 1 element
	/**
		 \relates MatrixOfElements
	 */		
	MatrixOfElements operator/(const MatrixOfElements & a, const VectorOfElements & b);

	
	///transposed matrix
	/**
		 \relates MatrixOfElements
	 */			
	MatrixOfElements transpose(MatrixOfElements & source);

	///element product of two vectors
	/**
		 \relates MatrixOfElements
		 \f$ elementProduct\left(
		                         \left[\begin{array}{c}
		                               a_1\\ \vdots \\ a_n
		                               \end{array}\right],
		                         \left[\begin{array}{c}
		                               b_1\\ \vdots \\ b_n
		                               \end{array}\right]\right) =
		                         \left[\begin{array}{ccc}
		                               a_1b_1 & \cdots & a_1b_n\\
		                               \vdots & \ddots & \vdots\\
		                               a_nb_1 & \cdots & a_nb_n\\
		                               \end{array}\right]
		                         
		 \f$, \f$A_{ij} = a_i b_j \f$ 
	 */			
	MatrixOfElements elementProduct(const VectorOfElements & a, const VectorOfElements & b);

	///Trace of a matrix \f$Tr(A)\equiv A_{ii}\f$ \relates MatrixOfElements		
	VectorOfElements trace(const MatrixOfElements & a); 

	///Trace of a matrix product \f$Tr(A B)\equiv A_{ij}B_{ji}\f$ \relates MatrixOfElements	
	VectorOfElements trace(const MatrixOfElements & a, const MatrixOfElements & b); 
	
	///	 generates a matrix with a row \relates MatrixOfElements
	MatrixOfElements generateME(const VectorOfElements & a);
	
	///	 generates a matrix with two rows 	 \relates MatrixOfElements
	MatrixOfElements generateME(const VectorOfElements & a,VectorOfElements & b);
	
	///	 generates a matrix with three rows  \relates MatrixOfElements
	MatrixOfElements generateME(const VectorOfElements & a,
	                            const VectorOfElements & b,
	                            const VectorOfElements & c);

	///	 generates a matrix with \p n rows \f$ generateME(\{u_i\}_j) = A_{ji}\f$ \relates MatrixOfElements
	MatrixOfElements generateME(const VectorOfElements *a, unsigned int n);

	///	 generates a matrix with \p n rows \f$ generateME(\{u_i\}_j) = A_{ji}\f$ \relates MatrixOfElements
	MatrixOfElements generateME(const vector<VectorOfElements> &a);
	
	/// returns VectorOfElements containing the diagonal elements
	/**
		\relates MatrixOfElements 
		the finction is valid only for square matrices
	*/
	VectorOfElements getDiagonal(const MatrixOfElements & a);
	
	/// returns VectorOfElements containing the uper off diagonal elements
	/**
		 \relates MatrixOfElements
		the finction is valid only for square matrices
	*/
	VectorOfElements getOffDiagonalUp(const MatrixOfElements & a);

	/// computes determinant expression fo cases 2x2 and 3x3 only
	/**
		 \relates MatrixOfElements		 
	*/
	VectorOfElements det(const MatrixOfElements & m);

	/// returns solution of a system of linear equations
	/**
		\ingroup ComplexDataTypes		 
	*/
	VectorOfElements solveSystem(const MatrixOfElements & a, const VectorOfElements & b);

	/// generates code for solving the solution of a system of linear equations
	/**
		\ingroup ComplexDataTypes
		Solves system of  2 or 3 equations with the Cramer's rule  
	*/
	vector<Element> gcSolveSystem(const MatrixOfElements & a, 
	                              const VectorOfElements & b, 
	                              const VectorOfElements & x);

	/// generates code for solving the solution of a system of linear equations
	/**
		\ingroup ComplexDataTypes		 
		 Solves system of  an arbiterary number of equations with the conjugate gradient method. 
		 The function constructs the algorithm with number of iterations of vector size. The initial
		 value of \f$ x_0 \f$ is \f$b\f$.
	*/
	vector<Element> gcSolveSystemCG(const MatrixOfElements & a, 
	                                const VectorOfElements & b, 
	                                const VectorOfElements & x);
	
	/// generate matrix with content of the matrix \p a but with replaced row \p r by vector \p b \relates MatrixOfElements		 
	MatrixOfElements replaceRow(const MatrixOfElements & a, const VectorOfElements & b, unsigned int r);

	/// generate matrix with content of the matrix \p a but with replaced column \p c by vector \p b \relates MatrixOfElements		 
	MatrixOfElements replaceColumn(const MatrixOfElements & a, const VectorOfElements & b, unsigned int c);

	/// returns the matrix of cofactors for cases 2x2 and 3x3 \relates MatrixOfElements 
	MatrixOfElements  generateMatrixCofactors(const MatrixOfElements & a);

	/// returns vector of elements for computing the inverse matrix for cases 2x2 and 3x3 \relates MatrixOfElements 
	vector<Element> gcMatrixInversion(const MatrixOfElements & a, MatrixOfElements & inv);

// ------------------------------ Implementation -----------------

	inline void MatrixOfElements::resize(unsigned int nr, unsigned int nc)
	{
		nRow=nr; 
		nCol=nc;
		ve.resize(nr*nc);
	}
	
	
}  //namespace acl

#endif // ACLMATRIXOFELEMENTS_H
