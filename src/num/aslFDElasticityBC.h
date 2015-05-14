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


#ifndef ASLFDELASTICITYBC_H
#define ASLFDELASTICITYBC_H

#include "aslBCond.h"
#include "acl/aclMath/aclVectorOfElementsDef.h"

namespace acl{
	class Kernel;
	typedef std::shared_ptr<Kernel> SPKernel;
}

namespace asl
{
	class ElasticityCommonA;
	typedef std::shared_ptr<ElasticityCommonA> SPElasticityCommonA;
	class FDElasticityIncompressibleStatic;
	typedef std::shared_ptr<FDElasticityIncompressibleStatic> SPFDElasticityIncompressibleStatic;
	class FDElasticityRelaxation;
	typedef std::shared_ptr<FDElasticityRelaxation> SPFDElasticityRelaxation;
	class FDElasticity2;
	typedef std::shared_ptr<FDElasticity2> SPFDElasticity2;
		
	/// Bondary condition corresponding to a rigid wall (\f$\vec u=0\f$ and \f$\nabla p=0\f$)
	/**	 
		 \ingroup ElasticityBC
	*/
	class BCRigidWall:public BCond
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;			
			SPFDElasticityIncompressibleStatic num;
		public:
			BCRigidWall(SPFDElasticityIncompressibleStatic nm);
			virtual void execute();
			virtual void init();			
	};

	/// Bondary condition corresponding to a rigid wall (\f$\vec u=0\f$ and \f$\nabla p=0\f$)
	/**	 
		 \ingroup ElasticityBS
	*/

	SPBCond generateBCRigidWall(SPFDElasticityIncompressibleStatic nm, 
                                 const std::vector<SlicesNames> & sl);

	/// Bondary condition corresponding to a rigid wall (\f$\vec u=0\f$ and \f$\nabla p=0\f$)
	/**	 
		 \ingroup ElasticityBC
	*/
	class BCRigidWallRelaxation:public BCond
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;			
			SPFDElasticityRelaxation num;
			acl::VectorOfElements value;
		public:
			BCRigidWallRelaxation(SPFDElasticityRelaxation nm);
			BCRigidWallRelaxation(SPFDElasticityRelaxation nm, acl::VectorOfElements v);
			virtual void execute();
			virtual void init();			
	};

	/// Bondary condition corresponding to a rigid wall (\f$\vec u=0\f$ and \f$\nabla p=0\f$)
	/**	 
		 \ingroup ElasticityBS
	*/

	SPBCond generateBCRigidWall(SPFDElasticityRelaxation nm, 
                                 const std::vector<SlicesNames> & sl);
	
	/// Bondary condition corresponding to a rigid wall (\f$ u=0\f$ and \f$\dot u=0\f$)
	/**	 
		 \ingroup ElasticityBC
	*/

	SPBCond generateBCRigidWall(SPFDElasticity2 nm, 
                                 const std::vector<SlicesNames> & sl);
	
	/// Bondary condition corresponding to a rigid wall (\f$ u=u_0\f$ and \f$\dot u=0\f$)
	/**	 
		 \ingroup ElasticityBC
	*/

	SPBCond generateBCRigidWall(SPFDElasticity2 nm, 
	                             const AVec<> & u0,
                                 const std::vector<SlicesNames> & sl);

	/// Bondary condition corresponding to a rigid wall (\f$ u=u_0\f$ and \f$\dot u=0\f$)
	/**	 
		 \ingroup ElasticityBC
	*/

	SPBCond generateBCRigidWall(SPFDElasticityRelaxation nm, 
	                             const AVec<> & u0,
                                 const std::vector<SlicesNames> & sl);
	
	/// Bondary condition corresponding to a free surface (\f$\partial_{\vec n}u=0\f$ and \f$\partial_{\vec n}\dot u=0\f$)
	/**	 
		 \ingroup ElasticityBC
	*/
	class BCFreeSurface:public BCond
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;			
			FDElasticityIncompressibleStatic* num;
		public:
			BCFreeSurface(FDElasticityIncompressibleStatic* nm);
			virtual void execute();
			virtual void init();			
	};

	/// Bondary condition corresponding to a free surface (\f$\partial_{\vec n}u=0\f$)
	/**	 
		 \note The free boundary conditions have a different mathematical formulation, 
		 see BCFreeSurface2Real 
		 
		 \ingroup ElasticityBC
	*/
	class BCFreeSurface2:public BCond
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;			
			SPFDElasticity2 num;
		public:
			BCFreeSurface2(SPFDElasticity2 nm);
			~BCFreeSurface2();			
			virtual void execute();
			virtual void init();			
	};

	/// Bondary condition corresponding to a free surface 
	/**	 
		 \ingroup ElasticityBC
		 
		The class is a realisation of 1st order explicit free boundary conditions
		for rectangular grid and arbitrary geometry. The geometry is defined by a levelset
		The algorithm computes value in a ghost node. This is done by cell subdivision. 
		The unit cell containing the ghost node in the center. It is splited into 
		elementary cells the spliting is described by the corresponding VTObjects. 

		 The gradient in the elementary cells is computed by the following expression:
		 \f[
		    \nabla_i u_j = (u^k_j - u^0_j E_k) T_{ki},
		 \f]
		 where \f$ u^k_j\f$ is the \f$ k^{th}\f$ point of the elementary cell, 
		 \f$ E_k \f$ is the vector filled by 1 and \f$ T_{ki} \f$ is the baricentric 
		 coordinates transformation matrix. 

		 The free boundary condition has the following expression:
		\f[
			\lambda n_i\nabla_k u_k+\mu n_j \left( \nabla_i u_j + \nabla_j u_i \right) = F_i,
		\f]
		where \f$ \lambda \f$ is a Lame first parameter, \f$ \mu \f$  is shear modulus.
		This expression can be rewriten in terms of the elementary cell expresion:
		\f[
			\lambda n_i(u^k_j - u^0_j E_k) T_{kj} + 
			\mu n_j \left( (u^k_j - u^0_j E_k) T_{ki} + (u^k_i - u^0_i E_k) T_{kj} \right) = F_i,
		\f]
		This equation can be rewritten in form of \f$ A_{ij} u^0_j = b_i \f$:
		\f[
		   A_{ij}=\lambda n_i E_kT_{kj} + \mu n_j  E_k T_{ki} + \mu n_l E_k T_{kl} \delta_{ij},
		\f]
		\f[
		   b_i = \lambda n_i u^k_j T_{kj} + \mu n_j (u^k_j T_{ki}+ u^k_i T_{kj}) - F_i.
		\f]		
		The solution of the system of equations is the seeking boundary condition.
			
		The class should obtain list of the ghost points

	*/
	class BCZeroStressMap: public BCondWithMap
	{		
		protected:
			SPAbstractDataWithGhostNodes displacement; 
			acl::VectorOfElements lambda;
			acl::VectorOfElements mu;
			acl::SPKernel kernel;
		public:
			BCZeroStressMap(SPAbstractDataWithGhostNodes displacement, 
			                acl::VectorOfElements lambda,
			                acl::VectorOfElements mu,
			                SPAbstractDataWithGhostNodes map,
			                const VectorTemplate *const t);
			~BCZeroStressMap();			
			virtual void execute();
			virtual void init();			
	};

	/// Bondary condition set given values to displacement/velocity
	/**	 
		 \ingroup ElasticityBC
	*/
	SPNumMethod generateBCZeroStress(SPElasticityCommonA nm, 
	                                 SPAbstractDataWithGhostNodes map);


	/// Bondary condition set given values to displacement/velocity
	/**	 
		 \ingroup ElasticityBC
	*/
//	SPNumMethod generateBCZeroStress(SPFDElasticityIncompressibleStatic nm, 
//	                                 SPAbstractDataWithGhostNodes map);

	/// Bondary condition set given values to displacement/velocity
	/**	 
		 \ingroup ElasticityBC
	*/
//	SPNumMethod generateBCZeroStress(SPFDElasticityRelaxation nm, 
//	                                 SPAbstractDataWithGhostNodes map);
	
	/// Bondary condition set given values to pressure
	/**	 
		 \ingroup ElasticityBC
	*/
	SPNumMethod generateBCZeroStressP(SPFDElasticityIncompressibleStatic nm, 
	                                  SPAbstractDataWithGhostNodes map);

	/// Bondary condition set given values to displacement/velocity
	/**	 
		 \ingroup ElasticityBC
	*/
	SPNumMethod generateBCZeroStressP(SPFDElasticityRelaxation nm, 
	                                  SPAbstractDataWithGhostNodes map);
	
	/// Bondary condition set given values to displacement/velocity
	/**	 
		 \ingroup ElasticityBC
	*/
	class BCImposedDisplacementVelocityValue:public BCond
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;			
			SPFDElasticityIncompressibleStatic num;
			AVec<double> displacement;
			bool bDisplacement;
			AVec<double> velocity;
			bool bVelocity;
			bool initialized;
		public:
			BCImposedDisplacementVelocityValue(SPFDElasticityIncompressibleStatic nm);
			virtual void execute();
			virtual void init();			
			void setDisplacement(AVec<double> d);
			void setVelocity(AVec<double> v);
	};

	/// Bondary condition local force(acceleration)
	/**	 
		 \ingroup ElasticityBC
	*/
	class BCAccelerationSource2: public BCond
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;			
			FDElasticity2* num;
			acl::VectorOfElements acceleration;
			bool initialized;
		public:
			BCAccelerationSource2(FDElasticity2* nm);
			virtual void execute();
			virtual void init();			
			void setAcceleration(AVec<double> a);
	};
		 
} //asl

#endif //ASLFDELASTICITYBC_H
