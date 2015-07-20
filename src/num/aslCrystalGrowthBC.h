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


#ifndef ASLCRYSTALGROWTHBC_H
#define ASLCRYSTALGROWTHBC_H

#include "aslBCond.h"
#include <data/aslDataWithGhostNodes.h>
#include <acl/aclMath/aclVectorOfElementsDef.h>


namespace acl
{
	class Kernel;
	typedef std::shared_ptr<Kernel> SPKernel;
}

namespace asl
{
	class PositionFunction;
	typedef std::shared_ptr<PositionFunction> SPPositionFunction;

	
	/** Boundary condition that makes gradient proportional to a surface concentration
		\ingroup BoundaryConditions */		 
	class BCLinearGrowthMap:public BCondWithMap
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernel;			
		protected:
			Data  data;
			acl::VectorOfElements cEq;
			acl::VectorOfElements beta;
		public:
			BCLinearGrowthMap(Data d, 
			                  const acl::VectorOfElements & cEq,
			                  const acl::VectorOfElements & beta,			                  
			                  Data map, 
			                  const VectorTemplate *const t);
			BCLinearGrowthMap(Data d, 
			                  const acl::VectorOfElements & cEq,
			                  const acl::VectorOfElements & beta,			                  
			                  Data map, 
			                  Data computationalDomain, 
			                  const VectorTemplate *const t);
			~BCLinearGrowthMap();
			virtual void execute();
			virtual void init();			
	};

	/// Boundary condition that makes gradient proportional to a surface concentration 
	/**	
		 \ingroup BoundaryConditions
		 The boundary condition corresponds to different system with surface reaction. 

		 Fist system is temperature reliase during the crystal growth process from melt.
		 There growth velocity defined by the following equation
		 \f[ \vec v = \vec n \beta (T -T^{eq}) \f]
		 where \f$ \beta \f$ is a kinetic coefficient and \f$ T^{eq} \f$ is an equilibrium temperature.
		 This condition can be expressed in terms of flux for puarlly diffusion equation
		 \f[ \vec j = D \vec \nabla T, \f]
		 \f[ \vec j = v L, \f]
		 where \f$ L \f$ is latent heat. These two equation allow to write an expression for temperature
		 on the boundary as follows:
		 \f[ \beta^* (T-T^{eq}) = \vec n \cdot \vec \nabla T, \f]
		 where \f$ \beta^* \equiv \frac{\beta L}{D} \f$ it corresponds to \p ::beta.
		 In the discrete case the two eequations for flux can not be combines so simply. 
		 Fist of all let us introduce two help values surface area per ghost node and 
		 flux per ghost point. There is set of template directions dirrectod towards computational domain 
		 \f$ \Omega\f$. The flux per ghost point can be expressed as follows 
		 \f[ j^{gh} = \sum_{i \in \Omega} w_i(T_i - T_0)  \f] 	
		 It is seen that this equation doen not contains the surface area explicetely.
		 \f[ j^{gh} = S^{gh} \beta^* (T-T^{eq}) \f]. Note, \f$ T \f$ is temperature on the surface
		 It can be expressed as alinear combunation of \f$ T_0 \f$ and \f$ T_i \f$ like this:
		 \f[ T = A T_0 + \sum_{i \in \Omega} B_i T_i \f]
		 Than the \f$ T_0 \f$ can be computed as follows:
		 \f[ T_0 = \frac{ \sum_{i\in \Omega} (w_i  - S^{gh}\beta^* B_i ) T_i + S^{gh}\beta^* T^{eq} }
                        {S^{gh}\beta^* A + \sum_{i\in \Omega} w_i} \f]

		 The coefficients \f$ A \f$ and \f$ B_i \f$ can be chosen as follows. 
		 \f[ T = \left(\sum_{i \in \Omega}\vec n \cdot\vec a_i\right)^{-1}
		         \sum_{i \in \Omega} (T_0(1-x_i) + T_i x_i) \vec n \cdot\vec a_i
                       
		 \f] 				
		 thus 
		 \f[ A = \left(\sum_{i \in \Omega}\vec n \cdot\vec a_i\right)^{-1}
		         \sum_{i \in \Omega} (1-x_i) \vec n \cdot\vec a_i
                       
		 \f] 				
		 \f[ B_i = \left(\sum_{i \in \Omega}\vec n \cdot\vec a_i\right)^{-1}
		           x_i \vec n \cdot\vec a_i
                       
		 \f] 		

		 \note The function is realized for 1 component only!

*/
		
	class BCLinearGrowthMap1:public BCondWithMap
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernel;			
		protected:
			Data  data;
			acl::VectorOfElements cEq;
			acl::VectorOfElements beta;
		public:
			BCLinearGrowthMap1(Data d, 
			                  const acl::VectorOfElements & cEq,
			                  const acl::VectorOfElements & beta,			                  
			                  Data map, 
			                  const VectorTemplate *const t);
			BCLinearGrowthMap1(Data d, 
			                  const acl::VectorOfElements & cEq,
			                  const acl::VectorOfElements & beta,			                  
			                  Data map, 
			                  Data computationalDomain, 
			                  const VectorTemplate *const t);
			~BCLinearGrowthMap1();
			virtual void execute();
			virtual void init();			
	};

	/** Boundary condition that makes gradient proportional to a surface concentration, second order
	 \ingroup BoundaryConditions */		 
	class BCLinearGrowthMap2:public BCondWithMap
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernelCN;			
			acl::SPKernel kernelGN;			
		protected:
			Data  data;
			acl::VectorOfElements cEq;
			acl::VectorOfElements beta;
		public:
			BCLinearGrowthMap2(Data d, 
			                  const acl::VectorOfElements & cEq,
			                  const acl::VectorOfElements & beta,			                  
			                  Data map, 
			                  const VectorTemplate *const t);
			BCLinearGrowthMap2(Data d, 
			                  const acl::VectorOfElements & cEq,
			                  const acl::VectorOfElements & beta,			                  
			                  Data map, 
			                  Data computationalDomain, 
			                  const VectorTemplate *const t);
			~BCLinearGrowthMap2();
			virtual void execute();
			virtual void init();			
	};

	
	
	SPNumMethod generateBCLinearGrowth(SPAbstractDataWithGhostNodes d, 
                                       double cEq,
                                       double beta,
                                       SPAbstractDataWithGhostNodes map,
                                       const VectorTemplate *const t);
	
	SPNumMethod generateBCLinearGrowth(SPAbstractDataWithGhostNodes d, 
	                                   double cEq,
                                       double beta,
	                                   SPAbstractDataWithGhostNodes map,
	                                   SPAbstractDataWithGhostNodes computationalDomain,
	                                   const VectorTemplate *const t);

	SPNumMethod generateBCLinearGrowth2(SPAbstractDataWithGhostNodes d, 
	                                    double cEq,
	                                    double beta,
	                                    SPAbstractDataWithGhostNodes map,
	                                    const VectorTemplate *const t);
	
	SPNumMethod generateBCLinearGrowth2(SPAbstractDataWithGhostNodes d, 
	                                    double cEq,
	                                    double beta,
	                                    SPAbstractDataWithGhostNodes map,
	                                    SPAbstractDataWithGhostNodes computationalDomain,
	                                    const VectorTemplate *const t);
	
} //asl

#endif //ASLCRYSTALGROWTHBC_H
