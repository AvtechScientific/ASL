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


#ifndef ASLUVALUE_H
#define ASLUVALUE_H

#include "aslTimeStamp.h"
#include <memory>

namespace asl
{

	/// Updatable value. This class stores value and its TimeStamp \ingroup LDI
	template<typename T> class UValue
	{
		public:
			std::shared_ptr<T> p;
			TimeStamp ts;
			inline explicit UValue():p(new T){}
			inline explicit UValue(const T & a):p(new T(a)){ts.update();}
			/// updates UValue::ts automaticaly
			inline const T & operator=(const T& a){ts.update(); return (*p = a);}
			inline const T & v() const {return *p;};
			inline T & v() {return *p;};
	};

	/// \related UValue
	template<typename T> inline bool operator==(const asl::UValue<T> & a,const T & v);
	/// \related UValue
	template<typename T> inline bool operator!=(const asl::UValue<T> & a,const T & v);

//-------------------------- Implementation --------------------------


	template<typename T> inline bool operator==(const asl::UValue<T> & a,const T & v)
	{
		return *(a.p) == v;
	}

	template<typename T> inline bool operator!=(const asl::UValue<T> & a,const T & v)
	{
		return *(a.p) != v;
	}
	
} //namespace asl


#endif // ASLUVALUE_H
