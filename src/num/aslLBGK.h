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


#ifndef ASLLBGK_H
#define ASLLBGK_H

#include "aslSingleKernelNM.h"
#include <acl/aclHardware.h>
#include <math/aslVectors.h>
#include "acl/aclMath/aclVectorOfElementsDef.h"
#include "acl/aclMath/aclMatrixOfElements.h"

namespace asl
{
	class Block;
	class VectorTemplate;
	template <typename V> class DataWithGhostNodes;
	typedef DataWithGhostNodes<acl::VectorOfElements> DataWithGhostNodesACL;
	typedef DataWithGhostNodes<acl::VectorOfElementsData> DataWithGhostNodesACLData;
	typedef std::shared_ptr<DataWithGhostNodesACL> SPDataWithGhostNodesACL;
	typedef std::shared_ptr<DataWithGhostNodesACLData> SPDataWithGhostNodesACLData;
	class AbstractDataWithGhostNodes;
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;


	/// returns VectorOfElements with values of rho  
	/**	
		 \relates LBGK

		 \f[ \rho=\sum_i w_i f_i\f]
		 where \f$ w_i \f$ is defined by \p vt 
	*/
	acl::VectorOfElements computeRho(acl::VectorOfElements f, const VectorTemplate* vt);
		
	/// returns VectorOfElements with values of velocity  
	/**	
		\relates LBGK
		 
		\f[ \vec m= \sum_i w_i\vec a_i f_i \f]
		where \f$ w_i \f$ is defined by \p vt 
	*/
	acl::VectorOfElements computeMomentum(acl::VectorOfElements f, const VectorTemplate* vt);

	/// generates Vector Of Elements with inverce components according to \p vt
	acl::VectorOfElements generateInverceVector(acl::VectorOfElements f, const VectorTemplate* vt);

	acl::MatrixOfElements generateLBGKMatrix(acl::VectorOfElements nu);
	acl::MatrixOfElements generateMRTMatrix(acl::VectorOfElements nu);
		
	acl::MatrixOfElements generateDifKinMatrix(acl::VectorOfElements nu);
	
	/// Numerical method for fluid flow
	/**
		 \ingroup TransportProcesses
		 \ingroup NumMethods		 
	*/
	class LBGK: public SingleKernelNM
	{
		public:
			typedef SPDataWithGhostNodesACL Data;
			typedef SPDataWithGhostNodesACLData DataD;
			typedef acl::VectorOfElements Param;
			
			const VectorTemplate* vectorTemplate;
		protected:
			acl::VectorOfElementsData fPool;
			Data f;
			DataD v;
			DataD rho;

			std::shared_ptr<AVec<int>> fShifts;
			std::shared_ptr<AVec<int>> fShiftsIncrement;

			std::vector<acl::SPKernel> copyKernels;

			Param viscosity;
			Param deltat;
			Param force;
			Param omega;

			bool flagComputeVelocity;
			bool flagComputeRho;
			bool flagCompressible;
			
			void createData(Block b, acl::CommandQueue queue, acl::TypeID type);
			void createCopyKernels();
			/// contains classical moving procedure
			virtual void preProcessing();
			virtual void init0();
			
		public:			
			LBGK();
			LBGK(DataD v, Param nu, const VectorTemplate* vT);
			LBGK(Block b, Param nu, const VectorTemplate* vT, 
			     bool compVel=true, bool compRho=true, 
			     acl::CommandQueue queue = acl::hardware.defaultQueue);
			void setViscosity(Param nu);
			double getViscosity(unsigned int i = 0);
			/// sets angular velocity for Coriolis term in noninertial reference frame
			void setOmega(Param w);
			void setVectorTemplate(const VectorTemplate* vT);
			inline const VectorTemplate* getVectorTemplate() const;


			inline Data getF();
			inline DataD getRho();
			inline DataD getVelocity();

			inline void setCompressible(bool flag = true);
			inline const bool & getCompressible() const;
	};

	typedef std::shared_ptr<LBGK> SPLBGK;	


	///contains different kernels for preprocessing and posprocessing of data used by LBGK
	class LBGKUtilities
	{
		private:
			typedef acl::VectorOfElements Param;
			
			SPLBGK num;
			acl::SPKernel kernelComputeV;
			acl::SPKernel kernelComputeRho;
			acl::SPKernel kernelComputeRhoV;
			acl::SPKernel kernelInitF;

			Param velocity;
			Param rho;
			
		public:
			LBGKUtilities(SPLBGK lbgk);
			void computeRho();
			void computeVelocity();
			void computeRhoVelocity();
			void initF(Param rho, Param vel);
			/// dencity is suposed to be 1
			void initF(Param vel);
	};

	typedef std::shared_ptr<LBGKUtilities> SPLBGKUtilities;

	class LBGKTurbulence: public LBGK
	{
		public:
			LBGKTurbulence(DataD v, Param nu, const VectorTemplate* vT);
			LBGKTurbulence(Block b, Param nu, const VectorTemplate* vT, 
			               bool compVel=true, bool compRho=true, 
			               acl::CommandQueue queue = acl::hardware.defaultQueue);
			virtual void init0();
	};

	typedef std::shared_ptr<LBGKTurbulence> SPLBGKTurbulence;	
	
// ------------------------- Implementation ------------------------


	inline LBGK::Data LBGK::getF()
	{
		return f;
	}

	inline LBGK::DataD LBGK::getRho()
	{
		return rho;
	}

	inline LBGK::DataD LBGK::getVelocity()
	{
		return v;
	}
	
	inline const VectorTemplate* LBGK::getVectorTemplate() const
	{
		return vectorTemplate;
	}
	
	inline void LBGK::setCompressible(bool flag)
	{
		flagCompressible=flag;
	}

	inline const bool & LBGK::getCompressible() const
	{
		return flagCompressible;
	}
		
} // asl
#endif // ASLLBGK_H
