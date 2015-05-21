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


#ifndef ASLLSNORMALGROWTH_H
#define ASLLSNORMALGROWTH_H

#include "aslLevelSetLinear.h"

namespace asl
{	
	/// Numerical method which computes evolution of an interface
	/**
		 \ingroup InterfaceTracking
		 \ingroup NumMethods
		 \ingroup LevelSet		 
	 */
	class LSNormalGrowth: public LevelSetLinear
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPDataWithGhostNodesACL DataGen;
			
		protected:
			DataGen superSaturation;
			
			virtual void initVelocityComputation();
		public:			

			LSNormalGrowth();
			/**
			 \param d is points position
			 \param c is super saturation field
			 */
			LSNormalGrowth(Data df, DataGen c);
			~LSNormalGrowth();
	};

	typedef std::shared_ptr<LSNormalGrowth> SPLSNormalGrowth;
	
} //asl

#endif //ASLLSNORMALGROWTH_H
