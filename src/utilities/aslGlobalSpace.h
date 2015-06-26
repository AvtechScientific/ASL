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


#ifndef ASLGLOBALSPACE_H
#define ASLGLOBALSPACE_H


namespace asl{
	
	//The class describes global discrete space
	/**
		It is used for correct positioning of the objects and 
		synchronization of processes in time
	*/
	class GlobalSpace {
		public:
			typedef long int Type;
			typedef AVec<long int> VType;
		private:
			const double dt;
			const double dx;
		public:
			inline GlobalSpace():dx(1.),dt(1.){};
			inline double discreteToGlobal(Type t){return dt*t;}
			inline Type globalToDiscrete(double t){return floor(t/dt);}
			inline Avec<double> discreteToGlobal(VType r){return dt*r;}
			inline VType globalToDiscrete(AVec<double> t){return floor(t/dt);}
	};

}  //namespace acl

#endif // ASLGLOBALSPACE_H
