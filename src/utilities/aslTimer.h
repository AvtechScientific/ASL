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

#include <chrono>
#include <thread>

namespace asl{
using namespace std::chrono;

	/// \ingroup Utilities
	class Timer {
		private:
			high_resolution_clock::time_point realTimeMark;
			duration<double> realTimeElapsed;

			clock_t processorTimeMark;
			double processorTimeElapsed;
		public:
			inline Timer() : processorTimeElapsed(0), realTimeElapsed(0) {};
			inline void start() {processorTimeMark = clock(); realTimeMark = high_resolution_clock::now();}
			inline void stop() {processorTimeElapsed += (double)(clock() - processorTimeMark); realTimeElapsed += duration_cast<duration<double>>(high_resolution_clock::now() - realTimeMark);}
			inline const double realTime() const{return realTimeElapsed.count();}
			inline const double processorTime() const{return (double)processorTimeElapsed/CLOCKS_PER_SEC;}
			inline const double processorLoad() const{return processorTime()/realTime();}
			inline void reset() {processorTimeElapsed = 0; realTimeElapsed = duration<double>(0);}
			/// Returns estimated duration of the current task based on its current \p completeness [0..1]
			inline const double estimatedDuration(double completeness) {return realTimeElapsed.count()/completeness;}
			/// Returns estimated time till finishing current task based on its current \p completeness [0..1]
			inline const double estimatedRemainder(double completeness) {return (1. - completeness) * estimatedDuration(completeness);}
	};

	/// Blocks execution of the calling thread for the time \p span (in milliseconds)
	inline void sleep(unsigned int span)
	{
		std::this_thread::sleep_for(std::chrono::milliseconds(span));
	}

}  //namespace acl

#endif // ASLTIMER_H
