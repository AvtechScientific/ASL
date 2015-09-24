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


/// \file aslVectorsUtilities.h Vectors and lattices

#ifndef _aslVectors
#define _aslVectors


#include "aslVectors.h"
#include <cmath>


namespace asl {


/*
  ///  The function makes possible an output to a stream
  template <int I,typename T> inline std::ostream& operator<<(std::ostream &f,const Vec<I,T> &a) {
    for (int i(0);i<I-1;++i)f<<a[i]<<' ';  return f<<a[I-1];
  }

  template <int I,typename T> inline const Vec<I,T> & operator+=(Vec<I,T> &a,const Vec<I,T> &b){
    for (int i(0);i<I;++i) a[i]+=b[i];  return a;
  }
  template <int I,typename T> inline const Vec<I,T> & operator-=(Vec<I,T> &a,const Vec<I,T> &b){
    for (int i(0);i<I;++i) a[i]-=b[i];  return a;
  }
  template <int I,typename T> inline const int operator== (const Vec<I,T> &a,const Vec<I,T> &b) {
    int c(1);  for (int i(0);i<I;++i) c&=(a[i]==b[i]);  return c;
  }
//  template <int I,typename T> inline const Vec<I,T> operator/ (const Vec<I,T> &b,const T &a) {return (1./a)*b;}

  ///Product of all Vector components, volume of ND dimentional box
  template <int I,typename T> inline const T bvol (const Vec<I,T> &a) {T p(1);  for (int i(0);i<I;++i) p*=a[i];  return p;}
  ///cheks whether the Vector components are nonegative
  template <int I,typename T,class IT> inline bool nonneg(const Vec<I,T,IT> &a) {
    int c(1);  for (int i(0);i<I;++i) c&=(a[i]>=0);  return c;
  }


  template <int I,typename T> inline const Vec<I,T> operator+ (const Vec<I,T> &a,const Vec<I,T> &b) {
    Vec<I,T> c;  for (int i(0);i<I;++i) c[i]=a[i]+b[i];  return c;
  };
  template <int I,typename T> inline const Vec<I,T> operator- (const Vec<I,T> &a,const Vec<I,T> &b) {
    Vec<I,T> c;  for (int i(0);i<I;++i) c[i]=a[i]-b[i];  return c;
  };
  template <int I,typename T> inline const Vec<I,T> operator* (const T &a,const Vec<I,T> &b) {
    Vec<I,T> c;  for (int i(0);i<I;++i) c[i]=a*b[i];  return c;
  };
  template <int I,typename T> inline const T operator* (const Vec<I,T> &a,const Vec<I,T> &b) {
    T c=0;  for (int i(0);i<I;++i) c+=a[i]*b[i];  return c;
  };

//  inline Vec<1,int> vceil(const Vec<1> &a) {return Vec<1,int>((int)ceil(a.x()));}
//  inline Vec<1,int> vfloor(const Vec<1> &a) {return Vec<1,int>((int)floor(a.x()));}
*/

  ///creates from two Vectors a Vector with maximal components
/*  template <typename T> inline Vec<3,T> cast_max(const Vec<3,T> &a,const Vec<3,T> &b) {
    return Vec<3,T>(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()));
    }
  ///creates from two Vectors a Vector with maximal components
  template <typename T> inline Vec<3,T> cast_min(const Vec<3,T> &a,const Vec<3,T> &b) {
    return Vec<3,T>(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()));
  }
  template <typename T> inline Vec<3,T> vabs(const Vec<3,T> &a) {
    return Vec<3,T>(fabs(a.x()),fabs(a.y()),fabs(a.z()));
  }*/

//  inline Vec<3,int> vround(const Vec<3> &a) {return Vec<3,int>((int)round(a.x()),(int)round(a.y()),(int)round(a.z()));}
//  inline Vec<3,int> vtrunc(const Vec<3> &a) {return Vec<3,int>((int)a.x(),(int)a.y(),(int)a.z());}


  ///returns minimum from the Vector components
//  template <typename T> inline const T min(const Vec<1,T>& a){return a.x();}
  ///returns maximum from the Vector components
//  template <typename T> inline const T max(const Vec<1,T>& a){return a.x();}

  ///Component production of Vector a on components of b
//  template <typename T> inline const Vec<3,T> cprod(const Vec<3,T> &a,const Vec<3,T> &b) {
//    return Vec<3,T>(a.x()*b.x(),a.y()*b.y(),a.z()*b.z());
//  }

/*  inline Vec<1,int> vec(int x){return Vec<1,int>(x);}
  inline Vec<2,int> vec(int x,int y){return Vec<2,int>(x,y);}
  inline Vec<3,int> vec(int x,int y,int z){return Vec<3,int>(x,y,z);}
  inline Vec<1> vec(lFl x){return Vec<1>(x);}
  inline Vec<2> vec(lFl x,lFl y){return Vec<2>(x,y);}
  inline Vec<3> vec(lFl x,lFl y,lFl z){return Vec<3>(x,y,z);}
*/
		
}

#endif
