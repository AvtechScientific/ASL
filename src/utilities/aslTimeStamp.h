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


#ifndef ASLTIMESTAMP_H
#define ASLTIMESTAMP_H


namespace asl{

	/// \ingroup LDI	
	class TimeStamp {
		public:
			typedef long unsigned int TSType;
		private:
			static TSType tsTime;
			TSType stamp;
		public:
			inline TimeStamp():stamp(0){};
			inline void update(){stamp=++tsTime; tsTime=(tsTime==0?1:tsTime);}
			inline bool needUpdate(const TimeStamp & a) const{
				return (stamp<a.stamp) || (stamp> tsTime);
			}
//			friend bool operator<(const TimeStamp &a, const TimeStamp &b);
//			friend bool operator>(const TimeStamp &a, const TimeStamp &b);
	};

/*	inline bool operator<(const TimeStamp &a, const TimeStamp &b)
	{
		return a.stamp<b.stamp;
	}

	inline bool operator>(const TimeStamp &a, const TimeStamp &b)
	{
		return a.stamp>b.stamp;
	}
*/

}  //namespace acl

#endif // ASLTIMESTAMP_H
