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


#ifndef ASLLBGKBC_H
#define ASLLBGKBC_H

#include "aslBCond.h"
#include "acl/aclMath/aclVectorOfElementsDef.h"

namespace acl{
	class Kernel;
	typedef std::shared_ptr<acl::Kernel> SPKernel;
	class KernelMerger; 
	typedef std::shared_ptr<acl::KernelMerger> SPKernelMerger;
}

namespace asl
{
	class LBGK;
	typedef std::shared_ptr<LBGK> SPLBGK;
	class PositionFunction;
	typedef std::shared_ptr<PositionFunction> SPPositionFunction;


	class BCLBGKCommon:public BCond
	{
		protected:
			SPLBGK num;
			std::vector<acl::SPKernel> kernels;
			acl::SPKernelMerger km;

			AVec<int> directionGroupsShifts;
			AVec<int> directionGroupsSizes;

			void sortDirections();
		public:
			BCLBGKCommon(SPLBGK nm);
			virtual void execute();			
	};
	
	/// Bondary condition corresponding to a rigid wall (\f$ \vec u=0\f$)
	/**	 
		 \ingroup TransportProcessesBC
		 The class realizes simple bounce back boundary conditions  
		 \f[ \vec v =0,\; \nabla P = 0 \f]
	*/
	class BCNoSlip: public BCLBGKCommon
	{		
		public:
			BCNoSlip(SPLBGK nm);
			virtual void init();			
	};

	/// Bondary condition corresponding an in- or outflow boundary conditions with a given pressure
	/**	 
		 \ingroup TransportProcessesBC
		 \f[ \vec v =0,\; P = Const \f]		 
	*/
	class BCConstantPressure: public BCLBGKCommon
	{		
		protected:
			acl::VectorOfElements pressure;
		public:
			BCConstantPressure(SPLBGK nm, const acl::VectorOfElements & p);
			virtual void init();			
	};

	/// Bondary condition corresponding wall with given velocity for uncompressible
	/**	 
		 \ingroup TransportProcessesBC
		 velocity value is valid for both tagential and normal components.
 		 \f[ \vec v =\vec v_0,\; \nabla P = 0 \f]		 

		 \f[
		    f^{ghost}_i = f^{bulk}_{\tilde i} + 2 \frac{\vec v \dot c_i}{c^2_s}
		 \f]

	*/
	class BCConstantVelocity: public BCLBGKCommon
	{	
		protected:
			acl::VectorOfElements velocity;
		public:
			BCConstantVelocity(SPLBGK nm, const acl::VectorOfElements & v);
			virtual void init();			
	};

	/// Bondary condition corresponding wall with given velocity for uncompressible
	/**	 
		 \ingroup TransportProcessesBC

		 velocity value is valid for both tagential and normal components.
 		 \f[ \vec v =\vec v_0,\; \nabla P = P_0 \f]		 

	*/
	class BCConstantPressureVelocity: public BCLBGKCommon
	{	
		protected:
			acl::VectorOfElements pressure;
			acl::VectorOfElements velocity;
		public:
			BCConstantPressureVelocity(SPLBGK nm,
			                           const acl::VectorOfElements & p,
			                           const acl::VectorOfElements & v);
			virtual void init();			
	};

	/// Bondary condition corresponding to a rigid wall (\f$ \vec u=0\f$)
	/**
		\f[ \vec v =0,\; \nabla P = 0 \f] 
	*/

	class BCNoSlipMap:public BCondWithMap
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;
			SPLBGK num;
		public:
			BCNoSlipMap(SPLBGK nm, SPAbstractDataWithGhostNodes map);
			~BCNoSlipMap();			
			virtual void execute();
			virtual void init();			
	};

	///
	/**
		 \ingroup TransportProcessesBC
	 */
	class BCVelocityMap:public BCondWithMap
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;
			SPLBGK num;
			SPPositionFunction velocity;
			
		public:
			BCVelocityMap(SPLBGK nm, 
			              SPPositionFunction v, 
			              SPAbstractDataWithGhostNodes map);
			BCVelocityMap(SPLBGK nm, 
			              SPPositionFunction v, 
			              SPAbstractDataWithGhostNodes map,
			              SPAbstractDataWithGhostNodes computationalDomain);
			~BCVelocityMap();			
			virtual void execute();
			virtual void init();			
	};

	/**
		\ingroup TransportProcessesBC
	 */	 
	class BCConstantPressureVelocityMap:public BCondWithMap
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;
			SPLBGK num;
			acl::VectorOfElements pressure;
			acl::VectorOfElements velocity;
		public:
			BCConstantPressureVelocityMap(SPLBGK nm, 
			          	          acl::VectorOfElements p, 
			                      SPAbstractDataWithGhostNodes map);
			BCConstantPressureVelocityMap(SPLBGK nm, 
			          	          acl::VectorOfElements p,
			                      acl::VectorOfElements v,
			                      SPAbstractDataWithGhostNodes map);
			~BCConstantPressureVelocityMap();			
			virtual void execute();
			virtual void init();			
	};

	/// Set outflux corresponding to transport limitation of the deposition rate
	/**
		 \ingroup TransportProcessesBC
		 
		 \f[ J_{dep}−J_{subl}=\frac{J_{dep}−P_0\alpha}{1+\alpha F, \f]
		 where \f$\alpha \equiv \sum_{i\in \Omega} w_i\f$,  \f$ F \f$ is the limiting factor

	*/
	class BCTransportLimitedDepositionMap:public BCondWithMap
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;
			SPLBGK num;
			acl::VectorOfElements p0;
			acl::VectorOfElements limitingFactor;
		public:
			BCTransportLimitedDepositionMap(SPLBGK nm, 
			                                acl::VectorOfElements p,
			                                acl::VectorOfElements lF,
			                                SPAbstractDataWithGhostNodes map);
			~BCTransportLimitedDepositionMap();			
			virtual void execute();
			virtual void init();			
	};

	/// Set outflux corresponding to kinetics and transport limitations of the deposition rate
	/**
		 \ingroup TransportProcessesBC
		 \f[ J_{dep}−J_{subl}=\beta\frac{\rho−\rho_0}{1+ F}, \f]
		 where \f$\alpha \equiv \sum_{i\in \Omega} w_i\f$,  \f$ F \f$ is the limiting factor

	*/
	class BCKineticsLimitedDepositionMap:public BCondWithMap
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;
			SPLBGK num;
			acl::VectorOfElements p0;
			acl::VectorOfElements limitingFactor;
			acl::VectorOfElements beta;
		public:
			BCKineticsLimitedDepositionMap(SPLBGK nm, 
			                                acl::VectorOfElements p,
			                                acl::VectorOfElements lF,
			                                acl::VectorOfElements b,
			                                SPAbstractDataWithGhostNodes map);
			~BCKineticsLimitedDepositionMap();			
			virtual void execute();
			virtual void init();			
	};
	
	/**
	 	 \ingroup TransportProcessesBC

	 */
	class ComputeSurfaceFluxMap:public BCondWithMap
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;
			SPLBGK num;
			SPDataWithGhostNodesACLData fluxField;
		public:
			ComputeSurfaceFluxMap(SPLBGK nm,
			                        SPDataWithGhostNodesACLData fF,
			                        SPAbstractDataWithGhostNodes map);
			~ComputeSurfaceFluxMap();			
			virtual void execute();
			virtual void init();			
	};

	/**
		 \ingroup TransportProcessesBC
	 */
	class ComputeSurfaceForceMap:public BCondWithMap
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;
			SPLBGK num;
			SPDataWithGhostNodesACLData forceField;
		public:
			ComputeSurfaceForceMap(SPLBGK nm,
			                        SPDataWithGhostNodesACLData fF,
			                        SPAbstractDataWithGhostNodes map);
			~ComputeSurfaceForceMap();			
			virtual void execute();
			virtual void init();			
	};

	///	\f[ \vec v =0,\; \nabla P = 0 \f] \ingroup TransportProcessesBC  
	SPBCond generateBCNoSlip(SPLBGK nm, const std::vector<SlicesNames> & sl);	
	///	\f[ \vec v =\vec v_0,\; \nabla P = 0 \f] \ingroup TransportProcessesBC 
	SPBCond generateBCConstantVelocity(SPLBGK nm, AVec<> v, const std::vector<SlicesNames> & sl);	
	///	\f[ \vec v =0,\; P = Const \f] \ingroup TransportProcessesBC 
	SPBCond generateBCConstantPressure(SPLBGK nm, double p, const std::vector<SlicesNames> & sl);	
	///	\f[ \vec v =\vec v_0,\; P = P_0 \f] \ingroup TransportProcessesBC
	SPBCond generateBCConstantPressureVelocity(SPLBGK nm, 
	                                   double p,
	                                   AVec<> v,
	                                   const std::vector<SlicesNames> & sl);	
	///	\f[ \vec v = 0,\; P = Const \f] \ingroup TransportProcessesBC 
	SPNumMethod generateBCConstantPressure(SPLBGK nm, double p, SPAbstractDataWithGhostNodes map);	
	///	\f[ \vec v =\vec v_0,\; P = P_0 \f] \ingroup TransportProcessesBC 
	SPNumMethod generateBCConstantPressureVelocity(SPLBGK nm, 
	                                               double p,
	                                               AVec<> v,            
	                                               SPAbstractDataWithGhostNodes map);	
	///	\f[ \vec v =0,\; \nabla P = 0 \f] \ingroup TransportProcessesBC 
	SPNumMethod generateBCNoSlip(SPLBGK nm, SPAbstractDataWithGhostNodes map);
	///	\f[ \vec v =0 \f] for velocity field \ingroup TransportProcessesBC 
	SPNumMethod generateBCNoSlipVel(SPLBGK nmU, SPAbstractDataWithGhostNodes map);
	///	\f[ \nabla P = 0 \f] \ingroup TransportProcessesBC 
	SPNumMethod generateBCNoSlipRho(SPLBGK nmU, SPAbstractDataWithGhostNodes map);
	///	\f[ \vec v =\vec v_0,\; \nabla P = 0 \f]  \ingroup TransportProcessesBC 
	SPNumMethod generateBCVelocity(SPLBGK nm, 
	                               SPPositionFunction v, 
	                               SPAbstractDataWithGhostNodes map);	
	///	\f[ \vec v =\vec v_0,\; \nabla P = 0 \f] \ingroup TransportProcessesBC 
	SPNumMethod generateBCVelocity(SPLBGK nm, 
	                               SPPositionFunction v, 
	                               SPAbstractDataWithGhostNodes map,
	                               SPAbstractDataWithGhostNodes computationalDomain);	
	///	\f[ \vec v =\vec v_0,\; \nabla P = 0 \f] \ingroup TransportProcessesBC 
	SPNumMethod generateBCVelocityVel(SPLBGK nm, 
	                                  SPPositionFunction v, 
	                                  SPAbstractDataWithGhostNodes map);
	///	\ingroup TransportProcessesBC 
	SPNumMethod generateBCTransportLimitedDeposition(SPLBGK nm, 
	                                                 double p0,
	                                                 double limitingFactor,            
	                                                 SPAbstractDataWithGhostNodes map);	
	///	\ingroup TransportProcessesBC 
	SPNumMethod generateBCKineticsLimitedDeposition(SPLBGK nm, 
	                                                double beta,
	                                                double p0,
	                                                double limitingFactor,            
	                                                SPAbstractDataWithGhostNodes map);	
	///	\ingroup TransportProcessesBC 
	SPNumMethod generateComputeSurfaceFlux(SPLBGK nm, 
	                                       SPDataWithGhostNodesACLData fF, 
	                                       SPAbstractDataWithGhostNodes map);	

	///	\ingroup TransportProcessesBC 
	SPNumMethod generateComputeSurfaceForce(SPLBGK nm, 
	                                       SPDataWithGhostNodesACLData fF, 
	                                       SPAbstractDataWithGhostNodes map);	
	
} //asl

#endif //ASLBGKBC_H
