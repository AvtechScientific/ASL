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


#ifndef ASLLSFACETEDGROWTH_H
#define ASLLSFACETEDGROWTH_H

#include "aslLevelSetLinear.h"
#include <math/aslVectorsDynamicLength.h>

namespace acl
{
	class ExpressionContainer;
}

namespace asl
{	

	/// describes crystalographyly specific date
	/**
	*/
	class CrystallographicParameters
	{
		public:
			vector<AVec<>> directions;
			vector<double> betaSt;
			vector<double> betaDisl;
			double betaRough;
		private:
			acl::VectorOfElements directionsACL;
			acl::VectorOfElements betaStACL;
			acl::VectorOfElements betaDislACL;	
		public:

			/// generates code which results the number corresponds to nearest direction
			void directionCode(acl::VectorOfElements normal,
				               acl::VectorOfElements direction,
			                   acl::VectorOfElements cosTheta,
				               acl::ExpressionContainer & k);
			CrystallographicParameters (const vector<AVec<>> & dir, 
				                       const vector<double> & bSt,
				                       const vector<double> & bDisl,
				                       const double bRough);
			CrystallographicParameters (const double bRough);
			CrystallographicParameters ();
			acl::VectorOfElements velocity(acl::VectorOfElements supersaturation, 
				                           acl::VectorOfElements dir,
			                               acl::VectorOfElements sinTheta);
			acl::VectorOfElements velocity(const acl::VectorOfElements & supersaturation, 
				                           const acl::VectorOfElements & dir,
			                               const acl::VectorOfElements & sinTheta,
			                               const acl::VectorOfElements & stepVelocityLimit);
			acl::VectorOfElements stepFactor(const acl::VectorOfElements & dir,
			                                 const acl::VectorOfElements & position);
		
			/// adds new facet information
			/**
				 \param normal should be a unit vector
			*/
			void addFacet(const AVec<> normal, double bSt, double bDisl);
			/// initialize acl data
			void init(acl::TypeID type);
	};

	
	/// Numerical method which computes evolution of an interface with a crystalographic kinetics
	/**
		 \ingroup InterfaceTracking
		 \ingroup NumMethods
		 \ingroup LevelSet		 
	 */
	class LSFacetedGrowth: public LevelSetLinear
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPDataWithGhostNodesACL DataGen;

			CrystallographicParameters crystallography;
		protected:
			DataGen superSaturation;
			
			virtual void initVelocityComputation();
		public:			

			LSFacetedGrowth();
			/**
			 \param d is points position
			 \param c is super saturation field
			 */
			LSFacetedGrowth(Data df, DataGen c);
			~LSFacetedGrowth();
	};

	typedef std::shared_ptr<LSFacetedGrowth> SPLSFacetedGrowth;
	
} //asl

#endif //ASLLSFACETEDGROWTH_H
