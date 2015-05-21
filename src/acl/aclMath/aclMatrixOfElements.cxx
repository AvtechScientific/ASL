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


#include "aclMatrixOfElements.h"
#include "aclVectorOfElements.h"
#include "aclGenerators.h"
#include "acl.h"
#include "../aslUtilities.h"
#include <algorithm>


using namespace std;
using asl::errorMessage;

namespace acl
{
	MatrixOfElements::MatrixOfElements(unsigned int nR, unsigned int nC):
		nRow(nR),
		nCol(nC),
		ve(nR*nC)
	{
	}

	unsigned int MatrixOfElements::ij2i(unsigned int i, unsigned int j) const
	{
		return i*nCol+j;
	}
	
	void MatrixOfElements::setElement(unsigned int i, unsigned int j, Element a)
	{
		ve.at(ij2i(i,j))=a;
	}

	void MatrixOfElements::setRow(unsigned int r, const VectorOfElements & b)
	{
		if(nCol != b.size())
			errorMessage("Error: MatrixOfElements::setRow: size of b does not match number of columns");
		if(nRow < r)
			errorMessage("Error: MatrixOfElements::setRow: r larger than number of rows");
		for(unsigned int i(0); i < b.size(); ++i)
			setElement(r,i,b[i]);
	}

	void MatrixOfElements::setColumn(unsigned int c, const VectorOfElements & b)
	{
		if(nRow != b.size())
			errorMessage("Error: MatrixOfElements::setColumn: size of b does not match number of rows");
		if(nCol < c)
			errorMessage("Error: MatrixOfElements::setColumn: c larger than number of columns");
		for(unsigned int i(0); i < b.size(); ++i)
			setElement(i,c,b[i]);
	}
		
	const Element MatrixOfElements::getElement(unsigned int i, unsigned int j) const
	{
		return ve.at(ij2i(i,j));
	}

	const VectorOfElements MatrixOfElements::getVE(unsigned int i, unsigned int j) const
	{
		return subVE(ve, ij2i(i,j));
	}

	const unsigned int MatrixOfElements::getNColumns() const
	{
		return nCol;
	}

	const unsigned int MatrixOfElements::getNRows() const
	{
		return nRow;
	}
	
	VectorOfElements & MatrixOfElements::getInternalVector()
	{
		return ve;
	}

	const VectorOfElements & MatrixOfElements::getInternalVector() const
	{
		return ve;
	}

	MatrixOfElements MatrixOfElements::operator= (const MatrixOfElements & m)
	{
		MatrixOfElements res(nRow, nCol);
		copy (ve = m.ve, res.ve);
		return res;
	}
		
		
	void copy(const MatrixOfElements & source, MatrixOfElements & destination)
	{
		if((source.getNRows() != destination.getNRows()) || 
		   (source.getNColumns() != destination.getNColumns()))
			errorMessage("copy: matricess have different sizes");
		copy(source.getInternalVector(),destination.getInternalVector());		
	}

	MatrixOfElements operator+(const MatrixOfElements & a, const MatrixOfElements & b)
	{
		if (a.getNRows() != b.getNRows() || a.getNColumns() != b.getNColumns())
			errorMessage("operator+ - the sizes of two MatrixOfElements are incompatible");

		MatrixOfElements m(a.getNRows(),a.getNColumns());
		copy(a.getInternalVector()+b.getInternalVector(),m.getInternalVector());
		return m;			
	}

	MatrixOfElements operator-(const MatrixOfElements & a, const MatrixOfElements & b)
	{
		if (a.getNRows() != b.getNRows() || a.getNColumns() != b.getNColumns())
			errorMessage("operator- - the sizes of two MatrixOfElements are incompatible");

		MatrixOfElements m(a.getNRows(),a.getNColumns());
		copy(a.getInternalVector()-b.getInternalVector(),m.getInternalVector());
		return m;
	}
		
	MatrixOfElements operator*(const MatrixOfElements & a, const MatrixOfElements & b)
	{
		if (a.getNColumns() != b.getNRows() && 
		    !(a.getNColumns()==1 && a.getNRows()==1) &&
		    !(b.getNColumns()==1 && b.getNRows()==1))
			errorMessage("operator* - the sizes of two MatrixOfElements are incompatible");

		MatrixOfElements m;
		if (a.getNColumns()==1 && a.getNRows()==1)
		{
			m.resize(b.getNRows(), b.getNColumns());
			copy(b.getInternalVector() * a.getInternalVector(), 
			     m.getInternalVector());
		} 
		else
		{
			if (b.getNColumns()==1 && b.getNRows()==1)
			{
				m.resize(a.getNRows(), a.getNColumns());
				copy(a.getInternalVector() * b.getInternalVector(), 
				     m.getInternalVector());
			}
			else
			{
				m.resize(a.getNRows(),b.getNColumns());
				for (unsigned int i(0); i < a.getNRows(); ++i)
					for (unsigned int j(0);j<b.getNColumns();++j)
					{
						using namespace elementOperators;
						Element r(a.getElement(i,0)*b.getElement(0,j));
						for (unsigned int k(1);k<b.getNRows();++k){
							r=r+a.getElement(i,k)*b.getElement(k,j);
						}
						m.setElement(i,j,r);
					}
			}
		}
		return m;			
	}

	VectorOfElements operator*(const VectorOfElements & a, const MatrixOfElements & b)
	{
		if (a.size() != b.getNRows())
			errorMessage("operator* - the sizes of VectorOfElements and MatrixOfElements are incompatible");
			
		VectorOfElements v(b.getNColumns());
		for (unsigned int j(0);j<b.getNColumns();++j)
		{
			using namespace elementOperators;
			Element r(a.at(0)*b.getElement(0,j));
			for (unsigned int k(1);k<b.getNRows();++k)
			{
				r=r+a.at(k)*b.getElement(k,j);
			}
			v[j] = r;
		}
		return v;
	}

	VectorOfElements operator*(const MatrixOfElements & a, const VectorOfElements & b)
	{
		if (a.getNColumns() != b.size())
			errorMessage("operator* - the sizes of MatrixOfElements and VectorOfElements are incompatible");

			
		VectorOfElements v(a.getNRows());
		for (unsigned int i(0); i < a.getNRows(); ++i){
			using namespace elementOperators;
			Element r(a.getElement(i,0)*b.at(0));
			for (unsigned int k(1);k<a.getNColumns();++k){
				r=r+a.getElement(i,k)*b.at(k);
			}
			v[i] = r;
		}
		return v;
	}

	MatrixOfElements operator/(const MatrixOfElements & a, const VectorOfElements & b)
	{
		if (b.size() !=1)
			errorMessage("operator/ - the sizes of Vector of elements is not 1");

			
		MatrixOfElements m(a.getNRows(),a.getNColumns());
		copy(a.getInternalVector()/b,m.getInternalVector());
		return m;
	}
		
	MatrixOfElements elementProduct(const VectorOfElements & a, const VectorOfElements & b){
		MatrixOfElements m(a.size(),b.size());
		using namespace elementOperators;
		for (unsigned int i(0); i < a.size(); ++i)
			for (unsigned int j(0);j<b.size();++j)
				m.setElement(i,j,a[i]*b[j]);
		return m;
	}

	VectorOfElements trace(const MatrixOfElements &a) 
	{
    	VectorOfElements c(a.getVE(0,0)); 
		for (unsigned int i(1); i < a.getNRows(); ++i) 
			copy(c + a.getVE(i, i), c); 
		return c;
	}

	VectorOfElements trace(const MatrixOfElements & a, const MatrixOfElements & b) 
	{
    	VectorOfElements c(generateVEConstant(0.));
		for (unsigned int i(0); i < a.getNRows(); ++i)
			for (unsigned int j(0); j < a.getNColumns(); ++j)
				copy(c + a.getVE(i,j) * b.getVE(j,i), c);
		return c;
	}

		
	MatrixOfElements generateME(const VectorOfElements & a)
	{
		MatrixOfElements m(1,a.size());
		copy(a,m.getInternalVector());
		return m;
	}
		
	MatrixOfElements generateME(const VectorOfElements & a, 
	                            const VectorOfElements & b)
	{
		if(a.size() != b.size())
			errorMessage("generateME: vectors have different sizes" );
		MatrixOfElements m(2,a.size());
		copy(cat(a,b),m.getInternalVector());
		return m;
	}
	
	MatrixOfElements generateME(const VectorOfElements & a,
	                            const VectorOfElements & b,
	                            const VectorOfElements & c)
	{
		if((a.size() != b.size()) || (a.size() != c.size()))
			errorMessage("generateME: vectors have different sizes" );
		MatrixOfElements m(3,a.size());
		copy(cat(a,b,c),m.getInternalVector());
		return m;	
	}
		
	MatrixOfElements generateME(const VectorOfElements *a, unsigned int n)
	{
		for(unsigned int i(1); i<n; ++i)
			if((a[0].size() != a[i].size()))
				errorMessage("generateME: some vectors have different sizes");

		MatrixOfElements m(n, a[0].size());
		for(unsigned int i(0); i<n; ++i)
			for(unsigned int j(0); j<a[i].size(); ++j)
				m.setElement(i,j,a[i][j]);
		return m;
	}

	MatrixOfElements generateME(const vector<VectorOfElements> & a)
	{
		return generateME(&a[0], a.size());
	}
		
	VectorOfElements getDiagonal(const MatrixOfElements & a)
	{
		if(a.getNColumns()!=a.getNRows())
			errorMessage("Error: getDiagonal: the matrix is not square one");

		unsigned int n(a.getNColumns());
		VectorOfElements v(n);

		for(unsigned int i(0); i<n; ++i)
			v[i]=a.getElement(i,i);

		return v;
	}

	VectorOfElements getOffDiagonalUp(const MatrixOfElements & a)
	{
		if(a.getNColumns()!=a.getNRows())
			errorMessage("Error: getOffDiagonalUp: the matrix is not square one");

		unsigned int n(a.getNColumns());
		VectorOfElements v((n*n-n)/2);

		unsigned int k(0);
		for(unsigned int i(1); i<n; ++i)
			for(unsigned int j(1); i<=i; ++j)
			{
				v[k]=a.getElement(i,j);
				++k;
			}
		return v;
	}

	/// computes determinant expression fo cases 2x2 and 3x3 only
	VectorOfElements det(const MatrixOfElements & m)
	{
		if(m.getNColumns()!=m.getNRows())
			errorMessage("Error: det: the matrix is not square one");
		if(m.getNColumns()>3)
			errorMessage("Error: det: the matrix size is larger than 3");

		VectorOfElements v(1);
		
		if(m.getNColumns()==1)
			copy(m.getInternalVector(),v);
		if(m.getNColumns()==2)
			copy(m.getVE(0,0) * m.getVE(1,1) - m.getVE(0,1) * m.getVE(1,0),
			     v);
//			copy(mad(m.getVE(0,0), m.getVE(1,1), - m.getVE(0,1) * m.getVE(1,0)),
//			     v);
		if(m.getNColumns()==3)
/*			copy(mad(m.getVE(0,0) * m.getVE(1,1),
			         m.getVE(2,2),
					 mad(m.getVE(1,0) * m.getVE(2,1), 
					     m.getVE(0,2),
						 mad(m.getVE(0,1) * m.getVE(1,2), 
							 m.getVE(2,0),
							 -mad(m.getVE(0,2) * m.getVE(1,1),
							     m.getVE(2,0),
							     mad(m.getVE(1,0) * m.getVE(0,1), 
							         m.getVE(2,2), 
							         m.getVE(2,1) * m.getVE(1,2) * m.getVE(0,0)))))),
			     v);*/
			copy(m.getVE(0,0) * m.getVE(1,1) * m.getVE(2,2) +
			     m.getVE(1,0) * m.getVE(2,1) * m.getVE(0,2) +
			     m.getVE(0,1) * m.getVE(1,2) * m.getVE(2,0) -
			     m.getVE(0,2) * m.getVE(1,1) * m.getVE(2,0) -
			     m.getVE(1,0) * m.getVE(0,1) * m.getVE(2,2) -
			     m.getVE(2,1) * m.getVE(1,2) * m.getVE(0,0),
			     v);

		return v;
	}

	VectorOfElements solveSystem(const MatrixOfElements & a, 
	                             const VectorOfElements & b)
	{
		if(a.getNColumns()!=a.getNRows())
			errorMessage("Error: solveSystem: the matrix is not square one");
		if(a.getNColumns()>3)
			errorMessage("Error: solveSystem: the matrix size is larger than 3");
		if(a.getNColumns()!=b.size())
			errorMessage("Error: solveSystem: size of b does not match the size of matrix");

		VectorOfElements v(b.size());

		auto d(det(a));
		VectorOfElements db(b.size());
		for(unsigned int i(0); i < b.size(); ++i)
			db[i] = (det(replaceColumn(a,b,i)))[0];
			
		return db/d;		
	}

	vector<Element> gcSolveSystem(const MatrixOfElements & a, 
	                              const VectorOfElements & b,
	                              const VectorOfElements & x)
	{
		if(a.getNColumns()!=a.getNRows())
			errorMessage("Error: solveSystem: the matrix is not square one");
		if(a.getNColumns()>3)
			errorMessage("Error: solveSystem: the matrix size is larger than 3");
		if(a.getNColumns()!=b.size())
			errorMessage("Error: solveSystem: size of b does not match the size of matrix");

		vector<Element> res(0);
		res.push_back((subVE(x,0) =  1./det(a))[0]);

		VectorOfElements db(b.size());
		for(unsigned int i(1); i < b.size(); ++i)
			res.push_back((subVE(x,i) = subVE(x,0) * det(replaceColumn(a,b,i)))[0]);
		res.push_back((subVE(x,0) *= det(replaceColumn(a,b,0)))[0]);
		return res;		
	}

	vector<Element> gcSolveSystemCG(const MatrixOfElements & a, 
	                                const VectorOfElements & b,
	                                const VectorOfElements & x)
	{
		if(a.getNColumns()!=a.getNRows())
			errorMessage("Error: solveSystem: the matrix is not square one");
		if(a.getNColumns()!=b.size())
			errorMessage("Error: solveSystem: size of b does not match the size of matrix");
		if(a.getNColumns()!=x.size())
			errorMessage("Error: solveSystem: size of x does not match the size of matrix");

		TypeID type(getElementType(x));
		auto p(generateVEPrivateVariable(b.size(),type));
		auto r(generateVEPrivateVariable(b.size(),type));
		auto a_p(generateVEPrivateVariable(b.size(),type));
		auto alpha(generateVEPrivateVariable(1,type));
		auto rr(generateVEPrivateVariable(1,type));
	
		vector<Element> res(0);
		res << (x=b);
		res << (r=b-a*x);
		res << (p=r);
		res << (rr = r*r);
		for(unsigned int i(0); i<b.size()-1; ++i)
		{
			res << (a_p = a * p);
			res << (alpha=rr/(p*a_p));
			res << (x+=alpha*p);
			res << (r-=alpha*a_p);
			// this line corresponds to beta
			res << (alpha=1./rr);
			res << (rr=r*r);
			res << (p = r+ rr*alpha * p);
		}
		res << (a_p = a * p);
		res << (alpha=rr/(p*a_p));
		res << (x+=alpha*p);
		
		return res;		
	}

		
	MatrixOfElements replaceRow(const MatrixOfElements & a,
	                            const VectorOfElements & b, 
	                            unsigned int r)
	{
		MatrixOfElements m(a.getNRows(),a.getNColumns());
		copy(a,m);
		m.setRow(r,b);
		return m;
	}

	MatrixOfElements replaceColumn(const MatrixOfElements & a, 
	                               const VectorOfElements & b, 
	                               unsigned int c)
	{
		MatrixOfElements m(a.getNRows(),a.getNColumns());
		copy(a,m);
		m.setColumn(c,b);
		return m;
	}

		
	VectorOfElements matrixCofactor(const MatrixOfElements & a, 
	                                unsigned int r, 
	                                unsigned int c)
	{
		MatrixOfElements m(a.getNRows()-1, a.getNColumns()-1);
		for(unsigned int i(0), im(0); im < m.getNRows(); ++i, ++im)
		{
			i += (i == r) ? 1 : 0;
			for(unsigned int j(0), jm(0); jm < m.getNColumns(); ++j, ++jm)
			{
				j += (j == c) ? 1 : 0;
				m.setElement(im, jm, a.getElement(i,j));
			}
		}
		return det(m);			
	}
		
	MatrixOfElements generateMatrixCofactors(const MatrixOfElements & a)
	{
		if(a.getNRows() != a.getNColumns())
			errorMessage("Error: generateMatrixCofactors: the matrix is not rectangular one");
		if(a.getNRows() > 3)
			errorMessage("Error: generateMatrixCofactors: the matrix size is more than 3");
		MatrixOfElements m(a.getNRows(), a.getNColumns());
		int rfactor(-1);
		for(unsigned int i(0); i < a.getNRows(); ++i){
			rfactor*=-1;
			int factor(-rfactor);
			for(unsigned int j(0); j < a.getNColumns(); ++j)
			{	
				factor*=-1;
				m.setElement(i,j, (factor * matrixCofactor(a,i,j))[0] );
			}
		}
		return m;
	}

	/// returns vector of elements for computing the inverse matrix for cases 2x2 and 3x3 \relates MatrixOfElements 
	vector<Element> gcMatrixInversion(const MatrixOfElements & a, MatrixOfElements & inv)
	{
		if((a.getNRows() != inv.getNRows()) || 
		   (a.getNColumns() != inv.getNColumns()))
			errorMessage("Error: generateMatrixInversionCode: two matrices have different sizes");
		auto mc(generateMatrixCofactors(a));
		auto d(det(a));

		unsigned int nr(a.getNRows());
		unsigned int nc(a.getNColumns());

		auto lastEl(inv.getVE(nr-1, nc-1));
		
		auto op1(lastEl = 1./det(a));
		auto op2(subVE(inv.getInternalVector(), 0, nr*nc-2) = 
		         subVE(mc.getInternalVector(), 0, nr*nc-2) * lastEl);
		auto op3(lastEl*=mc.getVE(nr-1, nc-1));

		return cat(op1,op2,op3);
	}
		
		
} // namespace acl
