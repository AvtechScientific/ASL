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


#ifndef ASLFDPOROELASTICITY_H
#define ASLFDPOROELASTICITY_H

#include "aslNumMethod.h"
#include "acl/aclMath/aclVectorOfElementsDef.h"
#include "aslFDElasticity.h"
#include "utilities/aslUValue.h"

namespace asl
{
	
	/// Numerical method which computes homogenious isotropic poro-elasticity equation
	/**
		 \ingroup Elasticity
		 \ingroup NumMethods

		The classic poroelastic equations were originally developed by Biot (1941) 
		to represent biphasic soil consolidation. The equations can be written 
		in the following form:
		\f[ (K+\mu/3)\nabla_j  \nabla_k u_k+ \mu \Delta u_j - a \nabla_j p = - (\rho_s-\rho_f)g_j  \f]
		\f[ a\nabla_i\dot u_i+ \frac{1}{S} \dot p = \nabla_i k \nabla_i p, \f]
		where 
			-  \f$\vec u\f$ is the displacement vector, 
			-  \f$ p \f$ the interstitial pressure,  
			-  \f$ a \f$ the ratio of fluid volume extracted to volume change of the tissue under compression, 
			-  \f$\rho_s\f$ the solid fraction density, 
			-  \f$\rho_f\f$ the fluid density, 
			-  \f$\vec g \f$ the gravitational vector, 
			-  \f$ 1/S \f$ the amount of fluid which can be forced into the tissue under constant volume, 
			-  \f$ k \f$ the hydraulic conductivity, 
			-  \f$K\f$ is the bulk modulus, 
			-  \f$\mu\f$ is the shear modulus,.

		In order to solve the first equation we will introduce a fictisionus time \f$ \tau \f$:
		\f[ \partial_\tau u = (K+\mu/3)\nabla_j  \nabla_k u_k+ \mu \Delta u_j - a \nabla_j p + (\rho_s-\rho_f)g_j  \f]
		\f[ a\partial_t\nabla_i u_i+ \frac{1}{S} \partial_t p = \nabla_i k \nabla_i p, \f]
			
			
		This implementation is aided to improove stability for incompressible materials.
		Let us reformulate the first eqation by the following way:
		\f[ \partial_\tau u = \mu \Delta u_j - \nabla_j \tilde p - a\nabla_j p + (\rho_s-\rho_f)g_j,  \f]
		\f[ \tilde p = -(K+\mu/3) \nabla_k u_k. \f]
		then the second eqation will turn to:
		\f[ \partial_t p = \frac{aS}{K+\mu/3}\partial_t\tilde p + S\nabla_i k \nabla_i p, \f]
			
		The equation for \f$\tilde p\f$ leads to an instability for low values of the Poisson ration (for nearly incompressible materials).
		Therefore the last equation we would like to replace by a relaxation one:
		\f[ \partial_\tau \tilde p = - w\left( (K+\mu/3)\nabla_k u_k + \tilde p\right), \f]
		where \f$ w \f$ is a relaxation parameter.

		Let's return to single time variable. Two time steps are related as follows 
			\f$ \delta_t = N \delta_\tau \f$. Than \f$ \partial_\tau = N\partial_t \f$. 
		The resulting set of equations is:
		\f[ \partial_\tau u = \mu \Delta u_j - \nabla_j \tilde p - a\nabla_j p + (\rho_s-\rho_f)g_j,  \f]
		\f[ \partial_\tau \tilde p = - w\left( (K+\mu/3)\nabla_k u_k + \tilde p\right), \f]			
		\f[ \partial_\tau p = \frac{aS}{K+\mu/3}\partial_\tau\tilde p + \frac{Sk}{N}\Delta p, \f]

	*/                    
	class FDPoroElasticity: public ElasticityCommonA
	{
		private:
			Data pressureData;
			Data pressureInternalData;			
			Data pressureLiquidData;
			Data pressureLiquidInternalData;

			Param hydraulicCondactivity;
			Param compresibility;
			Param nSubsteps;
				
		public:			
			FDPoroElasticity();
			/**
			 \param d is a displacement field
			 \param pl is a pressure of liquid field
			 \param bM is the bulk modulus
			 \param sM is the shear modulus
			 \param k is hydraulic conductivity 
			 \param vT is a vector template	
			 */
			FDPoroElasticity(Data d, Data pl, Param bM, 
			                 Param sM, Param k, 
			                 const VectorTemplate* vT);
			~FDPoroElasticity();

			virtual void init();
			virtual void execute();			

			inline Data getPressureData() const;
			inline Data getLiquidPressureData() const;

			/// defaul value 10
			void setNSubsteps(unsigned int n);
	};

	typedef std::shared_ptr<FDPoroElasticity> SPFDPoroElasticity;

	/**
	 \param d is a displacement field
	 \param pl is a pressure of liquid field
	 \param bM is the bulk modulus
	 \param sM is the shear modulus
	 \param k is hydraulic conductivity 
	 \param vT is a vector template	
	 */
	SPFDPoroElasticity generateFDPoroElasticity(SPDataWithGhostNodesACLData d,
	                                            SPDataWithGhostNodesACLData pl,
	                                            double bM, 
	                                            double sM,
	                                            double k,
	                                            const VectorTemplate* vT);
	
	/**
	 \param d is a displacement field
	 \param pl is a pressure of liquid field
	 \param bM is the bulk modulus
	 \param sM is the shear modulus
	 \param k is hydraulic conductivity 
	 \param vT is a vector template	
	 */
	template <typename T> 
		SPFDPoroElasticity generateFDPoroElasticity(SPDataWithGhostNodesACLData d,
		                                            SPDataWithGhostNodesACLData pl,
		                                            UValue<T> bM, 
		                                            UValue<T> sM,
		                                            UValue<T> k,
		                                            const VectorTemplate* vT);
	
//-------------------------IMPLEMENTATION------------------------

	inline ElasticityCommonA::Data FDPoroElasticity::getPressureData() const
	{
		return pressureData;
	}

	inline ElasticityCommonA::Data FDPoroElasticity::getLiquidPressureData() const
	{
		return pressureLiquidData;
	}
		
} // asl
#endif // ASLFDELASTICITY_H
