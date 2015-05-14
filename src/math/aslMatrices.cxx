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


#include"aslMatrices.h"

using namespace std;

namespace asl
{

	template <typename T> const AMatr<T> & AMatr<T>::operator=(const AMatr<T> & a)
	{
		if(nCol != a.nCol || nRow != a.nRow )
			resize(a.nRow, a.nCol);
		for(unsigned int i(0); i < nRow; ++i)
			for(unsigned int j(0); j < nCol; ++j)
				(*this)(i,j) = a(i,j);
		return *this;
	}
	
	template const AMatr<double> & AMatr<double>::operator=(const AMatr<double> & a);
	template const AMatr<float> & AMatr<float>::operator=(const AMatr<float> & a);
	
	template <typename T> void AMatr<T>::setRow(unsigned int r, const AVec<T> & b)
	{
		if(nCol != b.getSize())
			errorMessage("Error: AMatr<T>::setRow: size of b does not match number of columns");
		if(nRow < r)
			errorMessage("Error: AMatr<T>::setRow: r larger than number of rows");
		for(unsigned int i(0); i < b.getSize(); ++i)
			(*this)(r,i) = b[i];
	}

	template void AMatr<double>::setRow(unsigned int r, const AVec<double> & b);
	template void AMatr<float>::setRow(unsigned int r, const AVec<float> & b);
	
	template <typename T> void AMatr<T>::setColumn(unsigned int c, const AVec<T> & b)
	{
		if(nRow != b.getSize())
			errorMessage("Error: AMatr<T>::setColumn: size of b does not match number of rows");
		if(nCol < c)
			errorMessage("Error: AMatr<T>::setColumn: c larger than number of columns");
		for(unsigned int i(0); i < b.getSize(); ++i)
			(*this)(i,c) = b[i];
	}

	template void AMatr<double>::setColumn(unsigned int c, const AVec<double> & b);
	template void AMatr<float>::setColumn(unsigned int c, const AVec<float> & b);
	
	template <typename T> ostream& operator<<(ostream &f, const AMatr<T> & a) 
	{
		for (unsigned int i(0); i < a.getNRow(); ++i){
			for (unsigned int j(0); j < a.getNCol() - 1; ++j)
				f << a(i, j) << ' ';
			f << a(i, a.getNCol()-1) << endl;
		};
		return f;
	}

	template ostream& operator<<(ostream &f,const AMatr<float> & a);
	template ostream& operator<<(ostream &f,const AMatr<double> & a);
	
	template <typename T> const AMatr<T> operator* (const AMatr<T> &a, const AMatr<T> & b) 
	{
		AMatr<T> c(a.getNRow(), b.getNCol()); 
		T s(0);
		for (unsigned int i(0); i < a.getNRow(); ++i)
			for (unsigned int j(0); j < b.getNCol(); ++j)
			{
				for (unsigned int k(0); k<a.getNCol(); ++k) 
					s+=a(i, k) * b(k, j);
				c(i, j)=s; 
				s=0; 
			}
		return c;
	}

	template const AMatr<double> operator* (const AMatr<double> &a, const AMatr<double> & b); 
	template const AMatr<float> operator* (const AMatr<float> &a, const AMatr<float> & b);

	template <typename T> const AVec<T> operator* (const AMatr<T> &a, const AVec<T> & b) 
	{
		AVec<T> c(a.getNRow()); 
		T s(0);
		for (unsigned int i(0); i < a.getNRow(); ++i){
			for (unsigned int k(0); k<a.getNCol(); ++k) 
				s+=a(i, k) * b[k];
			c[i]=s; 
			s=0; 
		}
		return c;
	}

	template const AVec<double> operator* (const AMatr<double> &a, const AVec<double> & b); 
	template const AVec<float> operator* (const AMatr<float> &a, const AVec<float> & b);

	template <typename T> const AVec<T> operator* (const AVec<T> &a, const AMatr<T> & b) 
	{
		AVec<T> c(b.getNCol()); 
		T s(0);
		for (unsigned int j(0); j < b.getNCol(); ++j)
		{
			for (unsigned int k(0); k<a.getSize(); ++k) 
				s+=a[k] * b(k, j);
			c[j]=s; 
			s=0; 
		}
		return c;
	}

	template const AVec<double> operator* (const AVec<double> &a, const AMatr<double> & b); 
	template const AVec<float> operator* (const AVec<float> &a, const AMatr<float> & b);
	
	template <typename T> const AMatr<T> operator* (const AMatr<T> &a, const T & b) 
	{
		AMatr<T> c(a.getNRow(), a.getNCol()); 
		c.getInternalVec() = a.getInternalVec() * b;
		return c;
	}

	template const AMatr<double> operator* (const AMatr<double> &a, const double & b); 
	template const AMatr<float> operator* (const AMatr<float> &a, const float & b);

	template <typename T> const AMatr<T> operator* (const T &a, const AMatr<T> & b) 
	{
		AMatr<T> c(b.getNRow(), b.getNCol()); 
		c.getInternalVec() = b.getInternalVec() * a;
		return c;
	}

	template const AMatr<double> operator* (const double &a, const AMatr<double> & b); 
	template const AMatr<float> operator* (const float &a, const AMatr<float> & b);
	
	template <typename T> AMatr<T> elementProduct(const AVec<T> & a, const AVec<T> & b)
	{
		AMatr<T> m(a.getSize(),b.getSize());
		for (unsigned int i(0); i < a.getSize(); ++i)
			for (unsigned int j(0);j<b.getSize();++j)
				m(i,j)= a[i]*b[j];
		return m;
	}

	template AMatr<double> elementProduct(const AVec<double> & a, const AVec<double> & b);
	template AMatr<float> elementProduct(const AVec<float> & a, const AVec<float> & b);
	
	template <typename T> const T trace(const AMatr<T> &a) 
	{
    	T c(0); 
		for (unsigned int i(0); i < a.getNRow(); ++i) 
			c += a(i, i); 
		return c;
	}
	
	template const double trace(const AMatr<double> &a); 
	template const float trace(const AMatr<float> &a); 

 	template <typename T> 
		const T trace(const AMatr<T> & a, const AMatr<T> & b) 
	{
    	T c(0);
		for (unsigned int i(0); i < a.getNRow(); ++i)
			for (unsigned int j(0); j < a.getNCol(); ++j)
				c += a(i,j) * b(j,i);
		return c;
	}

 	template const double trace(const AMatr<double> & a, const AMatr<double> & b); 
 	template const float trace(const AMatr<float> & a, const AMatr<float> & b); 

	template <typename T> AMatr<T> makeAMatr(const AVec<T> & a)
	{
		return {1,a.getSize(),a};
	}

	template AMatr<double> makeAMatr(const AVec<double> & a);
	template AMatr<float> makeAMatr(const AVec<float> & a);
		                       
	template <typename T> AMatr<T> makeAMatr(const AVec<T> & a, const AVec<T> & b)
	{
		if(a.getSize() != b.getSize())
			errorMessage("makeMatr: vectors have different sizes" );
		AMatr<T> m(2, a.getSize());
		m.setRow(0,a);
		m.setRow(1,b);		
		return m;
	}

	template AMatr<double> makeAMatr<double>(const AVec<double> & a, const AVec<double> & b);
	template AMatr<float> makeAMatr<float>(const AVec<float> & a, const AVec<float> & b);

	template <typename T> AMatr<T> makeAMatr(const AVec<T> & a, 
	                                         const AVec<T> & b, 
	                                         const AVec<T> & c)
	{
		if((a.getSize() != b.getSize()) || (a.getSize() != c.getSize()))
			errorMessage("makeMatr: vectors have different sizes" );
		AMatr<T> m(3, a.getSize());
		m.setRow(0,a);
		m.setRow(1,b);		
		m.setRow(2,c);		
		return m;
	}

	template AMatr<double> makeAMatr(const AVec<double> & a, const AVec<double> & b, const AVec<double> & c);
	template AMatr<float> makeAMatr(const AVec<float> & a, const AVec<float> & b, const AVec<float> & c);
			
	template <typename T> AMatr<T> makeAMatr(AVec<T> *a, unsigned int n)
	{
		for(unsigned int i(1); i<n; ++i)
			if((a[0].getSize() != a[i].getSize()))
				errorMessage("makeMatr: some vectors have different sizes");

		AMatr<T> m(n, a[0].getSize());
		for(unsigned int i(0); i<n; ++i)
				m.setRow(i,a[i]);
		return m;
	}

	template AMatr<double> makeAMatr(AVec<double> *a, unsigned int n);
	template AMatr<float> makeAMatr(AVec<float> *a, unsigned int n);

	template <typename T> AMatr<T> makeAMatrUnit(unsigned int n)
	{
		AMatr<T> m(n,n);
		
		for(unsigned int i(0); i<n; ++i)
			for(unsigned int j(0); j<n; ++j)
				m(i,j)=0;
		for(unsigned int i(0); i<n; ++i)
			m(i,i) = 1;
		return m;
	}

	template  AMatr<double> makeAMatrUnit<double>(unsigned int n);
	template  AMatr<float> makeAMatrUnit<float>(unsigned int n);
	
	template <typename T> AVec<T> getDiagonal(const AMatr<T> & a)
	{
		if(a.getNCol()!=a.getNRow())
			errorMessage("Error: getDiagonal: the matrix is not square one");

		unsigned int n(a.getNCol());
		AVec<T> v(n);

		for(unsigned int i(0); i<n; ++i)
			v[i]=a(i,i);

		return v;
	}

	template AVec<double> getDiagonal(const AMatr<double> & a);
	template AVec<float> getDiagonal(const AMatr<float> & a);
	
	template <typename T> AVec<T> getOffDiagonalUp(const AMatr<T> & a)
	{
		if(a.getNCol()!=a.getNRow())
			errorMessage("Error: getOffDiagonalUp: the matrix is not square one");

		unsigned int n(a.getNCol());
		AVec<T> v((n*n-n)/2);

		unsigned int k(0);
		for(unsigned int i(1); i<n; ++i)
			for(unsigned int j(1); i<=i; ++j)
			{
				v[k]=a.getElement(i,j);
				++k;
			}
		return v;
	}

	template <typename T> T det(const AMatr<T> & m)
	{
		if(m.getNCol()!=m.getNRow())
			errorMessage("Error: det: the matrix is not square one");
		if(m.getNCol()>3)
			errorMessage("Error: det: the matrix size is larger than 3");

		T v(0);
		if(m.getNCol()==1)
			v=m(0,0);
		if(m.getNCol()==2)
			v=m(0,0) * m(1,1) - m(0,1) * m(1,0);
		if(m.getNCol()==3)
			v=m(0,0) * m(1,1) * m(2,2) +
			  m(1,0) * m(2,1) * m(0,2) +
			  m(0,1) * m(1,2) * m(2,0) -
			  m(0,2) * m(1,1) * m(2,0) -
			  m(1,0) * m(0,1) * m(2,2) -
			  m(2,1) * m(1,2) * m(0,0);
		return v;
	}

	template double det(const AMatr<double> & m);
	template float det(const AMatr<float> & m);

	template <typename T> 
		AVec<T> solveSystem(const AMatr<T> & a, const AVec<T> & b)
	{
		if(a.getNCol()!=a.getNRow())
			errorMessage("Error: solveSystem: the matrix is not square one");
		if(a.getNCol()>3)
			errorMessage("Error: solveSystem: the matrix size is larger than 3");
		if(a.getNCol()!=b.getSize())
			errorMessage("Error: solveSystem: size of b does not match the size of matrix");

		AVec<T> v(b.getSize());

		T d(det(a));
		AVec<T> db(b.getSize());
		for(unsigned int i(0); i < b.getSize(); ++i)
			db[i] = det(replaceColumn(a,b,i));
			
		return db/d;		
	}

	template AVec<double> solveSystem(const AMatr<double> & a, const AVec<double> & b);
	template AVec<float> solveSystem(const AMatr<float> & a, const AVec<float> & b);
	
	
	template <typename T> AMatr<T> replaceRow(const AMatr<T> & a,
	                            const AVec<T> & b, 
	                            unsigned int r)
	{
		AMatr<T> m(a);
		m.setRow(r,b);
		return m;
	}

	template <typename T> AMatr<T> replaceColumn(const AMatr<T> & a, 
	                               const AVec<T> & b, 
	                               unsigned int c)
	{
		AMatr<T> m(a);
		m.setColumn(c,b);
		return m;
	}

		
	template <typename T> 
		T matrixCofactor(const AMatr<T> & a, unsigned int r, unsigned int c)
	{
		AMatr<T> m(a.getNRow()-1, a.getNCol()-1);
		for(unsigned int i(0), im(0); im < m.getNRow(); ++i, ++im)
		{
			i += (i == r) ? 1 : 0;
			for(unsigned int j(0), jm(0); jm < m.getNCol(); ++j, ++jm)
			{
				j += (j == c) ? 1 : 0;
				m(im, jm) = a(i,j);
			}
		}
		return det(m);			
	}

	
	template <typename T> AMatr<T> generateMatrixCofactors(const AMatr<T> & a)
	{
		if(a.getNRow() != a.getNCol())
			errorMessage("Error: generateMatrixCofactors: the matrix is not rectangular one");
		if(a.getNRow() > 3)
			errorMessage("Error: generateMatrixCofactors: the matrix size is more than 3");
		AMatr<T> m(a.getNRow(), a.getNCol());
		int rfactor(-1);
		for(unsigned int i(0); i < a.getNRow(); ++i){
			rfactor*=-1;
			int factor(-rfactor);
			for(unsigned int j(0); j < a.getNCol(); ++j)
			{	
				factor*=-1;
				m(i,j) = (factor * matrixCofactor(a,i,j));
			}
		}
		return m;
	}

	template <typename T> AMatr<T> inverseMatrix(const AMatr<T> & a)
	{
		AMatr<T> mc(generateMatrixCofactors(a));
		T d(det(a));
		return mc/d;
	}
	
	template AMatr<double> inverseMatrix(const AMatr<double> & a);
	template AMatr<float> inverseMatrix(const AMatr<float> & a);

	
	
	double detSymMat(double a,double b, double c)
	{
		return a*b-c*c;
	}

	void getEValEVecMatSym2x2(double a,double b, double c, double & l1, double & l2, double & v1x, double & v1y, double & v2x, double & v2y)
	{
		double tr2((a+b)*.5);
		double det(detSymMat(a,b,c));
		l1=	tr2 + sqrt(tr2*tr2-det);
		l2=	tr2 - sqrt(tr2*tr2-det);
		v1x=l1-b;
		v1y=c;
		v2x=l2-b;
		v2y=c;
	}

	double detSymMat(double a, double b, double c, double d, double e, double f)
	{
		return a*b*c+2.*d*e*f-b*e*e-a*f*f-c*d*d;
	}

	void getEValEVecMatSym3x3(double a, double b, double c, double d, double e, double f, 
		                      double & l1, double & l2, double & l3, 
		                      double & v1x, double & v1y, double & v1z, 
		                      double & v2x, double & v2y, double & v2z,
		                      double & v3x, double & v3y, double & v3z)
	{
//	double tr(a+b+c);	
//	double det(detSymMat(a,b,c,d,e,f))
		
	}

}
