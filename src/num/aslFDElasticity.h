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


#ifndef ASLFDELASTICITY_H
#define ASLFDELASTICITY_H

#include "aslNumMethod.h"
#include "acl/aclMath/aclVectorOfElementsDef.h"

namespace acl{
	class Kernel;
	class VectorOfElementsData;
}

namespace asl
{
	class VectorTemplate;
	class TemplateVE;

	template <typename V> class DataWithGhostNodes;
	typedef DataWithGhostNodes<acl::VectorOfElementsData> DataWithGhostNodesACLData;
	typedef std::shared_ptr<DataWithGhostNodesACLData> SPDataWithGhostNodesACLData;
	typedef DataWithGhostNodes<acl::VectorOfElements> DataWithGhostNodesACL;
	typedef std::shared_ptr<DataWithGhostNodesACL> SPDataWithGhostNodesACL;

	namespace elasticity
	{
		acl::VectorOfElements strain(vector<TemplateVE> & displacment);
		acl::VectorOfElements stressLinear(acl::VectorOfElements & strain);
//		VectorOfElements linearEqLaplasTerm(TemplateVE displacment);
	}

	/// abstract class for elasticity solver
	class ElasticityCommonA: public NumMethod
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef acl::VectorOfElements Param;
		protected:
			std::unique_ptr<acl::Kernel> kernel;

			Data displacementData;
			Data displacementInternalData;			
				
			Param bulkModulus;
			Param shearModulus;
			Param force;
		public:			
			const VectorTemplate* vectorTemplate;

			ElasticityCommonA();
			/**
			 \param d is a displacement field
			 \param bM is the bulk modulus
			 \param sM is the shear modulus
			 \param vT is a vector template			 
			 */
			ElasticityCommonA(Data d, Param bM, Param sM, const VectorTemplate* vT);
			~ElasticityCommonA();
			void setVectorTemplate(const VectorTemplate* vT);
			VectorTemplate getVectorTemplate();

			virtual void init()=0;
			virtual void execute()=0;			

			void setForce(Param f);
			inline Data getDisplacementData() const;
			inline Data getDisplacementInternalData() const;
			inline const Param getBulkModulus() const;
			inline const Param getShearModulus() const;
	};
	
	/// Numerical method which computes homogenious isotropic elasticity equation
	/**
		 \ingroup Elasticity
		 \ingroup NumMethods
	 */
	class FDElasticityIncompressibleStatic: public ElasticityCommonA
	{
		private:
			Data pressure;
			Data pressureInternalData;
		public:			
			FDElasticityIncompressibleStatic();
			/**
			 \param d is a displacement field
			 \param bM is the bulk modulus
			 \param sM is the shear modulus
			 \param vT is a vector template			 
			 */
			FDElasticityIncompressibleStatic(Data d, Param bM, Param sM, const VectorTemplate* vT);
			~FDElasticityIncompressibleStatic();

			virtual void init();
			virtual void execute();			

			inline Data getPressureData() const;
	};

	typedef std::shared_ptr<FDElasticityIncompressibleStatic> SPFDElasticityIncompressibleStatic;

	/// Numerical method which computes homogenious isotropic elasticity equation
	/**
		 \ingroup Elasticity
		 \ingroup NumMethods
	 */
	class FDElasticityRelaxation: public ElasticityCommonA
	{
		private:
			Data pressure;
			Data pressureInternalData;
			Param deltat;
			Param dumpingFactor;
		public:			
			FDElasticityRelaxation();
			/**
			 \param d is a displacement field
			 \param bM is the bulk modulus
			 \param sM is the shear modulus
			 \param dt is time step
			 \param vT is a vector template			 
			 */
			FDElasticityRelaxation(Data d, Param bM, Param sM, Param dt, const VectorTemplate* vT);
			~FDElasticityRelaxation();

			virtual void init();
			virtual void execute();			

			inline Data getPressureData() const;
			inline Param getDeltat() const;
			void setDumpingFactor(Param dumpF);
	};

	typedef std::shared_ptr<FDElasticityRelaxation> SPFDElasticityRelaxation;

	
	/// Numerical method which computes homogenious isotropic elasticity equation
	/**
		 \ingroup Elasticity
		 \ingroup NumMethods
		 
		 \f[ \rho\ddot u_j =(K+\mu/3)\nabla_j  \nabla_k u_k+ \mu \Delta u_j + \vec F \f]
		 where \f$K\f$ is the bulk modulus, \f$\mu\f$ is the shear modulus, 
		 \f$\vec u\f$ is a displacement vector field, \f$ \vec F\f$ external force
	*/
	class FDElasticity2: public ElasticityCommonA
	{
		private:
//			VData velocityData;
			Param deltat;
			Param dumpingFactor;
		public:			
			FDElasticity2();
			/**
			 \param d is a displacement field
			 \param bM is the bulk modulus
			 \param sM is the shear modulus
			 \param dt is time step
			 \param vT is a vector template			 
			 */
			FDElasticity2(Data d, Param bM, Param sM, 
			              Param dt, const VectorTemplate* vT);
			~FDElasticity2();

			virtual void init();
			virtual void execute();			
//			inline VData getVelocityData() const;

			inline Param getDeltat() const;
			void setDumpingFactor(Param dumpF);
	};

	typedef std::shared_ptr<FDElasticity2> SPFDElasticity2;
	
	/// Number of dimensions
	/**
		 \relates FDElasticity
	 */
	unsigned int nD(const ElasticityCommonA & e);

	SPFDElasticity2 generateFDElasticity(SPDataWithGhostNodesACLData d, 
	                                     double bM, 
	                                     double sM, 
	                                     double dt, 
	                                     const VectorTemplate* vT);

	SPFDElasticityIncompressibleStatic generateFDElasticityStatic(SPDataWithGhostNodesACLData d, 
	                                                              double bM, 
	                                                              double sM, 
	                                                              const VectorTemplate* vT);

	SPFDElasticityRelaxation generateFDElasticityRelax(SPDataWithGhostNodesACLData d, 
	                                                   double bM, 
	                                                   double sM, 
	                                                   double dt, 
	                                                   const VectorTemplate* vT);
	
//-------------------------IMPLEMENTATION------------------------

	inline ElasticityCommonA::Data 
		ElasticityCommonA::getDisplacementData() const
	{
		return displacementData;
	}

	inline ElasticityCommonA::Data ElasticityCommonA::getDisplacementInternalData() const
	{
		return displacementInternalData;
	}
	
	inline const ElasticityCommonA::Param 
		ElasticityCommonA::getBulkModulus() const
	{
		return bulkModulus;
	}
	
	inline const ElasticityCommonA::Param 
		ElasticityCommonA::getShearModulus() const
	{
		return shearModulus;
	}

	inline ElasticityCommonA::Data 
		FDElasticityIncompressibleStatic::getPressureData() const
	{
		return pressure;
	}
	
	inline FDElasticityRelaxation::Data 
		FDElasticityRelaxation::getPressureData() const
	{
		return pressure;
	}
	
	inline FDElasticity2::Param FDElasticity2::getDeltat() const
	{
		return deltat;
	}
	
} // asl
#endif // ASLFDELASTICITY_H
