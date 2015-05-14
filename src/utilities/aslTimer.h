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


#ifndef ASLTIMER_H
#define ASLTIMER_H

#include <sys/time.h>

namespace asl{

	/// \ingroup Utilities
	class Timer {
		private:
			clock_t _c;
			double _t;
			inline double tod(){
				timeval tim;
				gettimeofday(&tim,NULL);
				return tim.tv_sec+(tim.tv_usec/1000000.0);
			}
		public:
			inline Timer():_c(0),_t(0){};
			inline void start(){_c=-clock();_t=-tod();}
			inline void resume(){_c-=clock();_t-=tod();}
			inline void stop(){_c+=clock();_t+=tod();}
			inline const double getTime() const{return _t;}
			inline const double getClockTime() const{return (float)_c/CLOCKS_PER_SEC;}
			inline const double getProcessorLoad() const{return getClockTime()/getTime();}
			inline void reset(){_c=0; _t=0;}
			inline const double getExpectedTime(double fractTime){return _t/fractTime;}
			inline const double getLeftTime(double fractTime){return (1.-fractTime)*getExpectedTime(fractTime);}
	};

	// waits \p dt seconds, uses loop
	inline void delay(double dt)
	{
		timeval t;
		gettimeofday(&t,NULL);
		double t0(t.tv_sec+(t.tv_usec/1000000.0));
		double tc(t0);
		while (t0+dt>tc)
		{
			gettimeofday(&t,NULL);
			tc=t.tv_sec+(t.tv_usec/1000000.0);			
		}
	 }

}  //namespace acl

#endif // aslGeneratorS_H
