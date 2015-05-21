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


#ifndef ASLFDSTEFANMAXWELL_H
#define ASLFDSTEFANMAXWELL_H

#include "aslSingleKernelNM.h"

namespace acl
{
	class VectorOfElementsData;
	class VectorOfElements;
}

namespace asl
{
	class VectorTemplate;
	template <typename V> class DataWithGhostNodes;
	typedef DataWithGhostNodes<acl::VectorOfElementsData> DataWithGhostNodesACLData;
	typedef std::shared_ptr<DataWithGhostNodesACLData> SPDataWithGhostNodesACLData;
	class AbstractDataWithGhostNodes;
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;
	
	/// Numerical method which computes multicomponent transport processes
	/**
		 \ingroup TransportProcesses
		 \ingroup NumMethods
		 
		 \f[ \partial_t c_i= - \vec \nabla \cdot \vec J - \vec \nabla \cdot (\vec v c_i) \f]
		 \f[ -\nabla c_i = \sum_{j, i\neq j} \frac{c_j\vec J_i-c_i\vec J_j}{c_tD_{ij}} + \frac{\vec J_i}{D_{i,D}}\f]						
		 where \f$c_i\f$ is a molar concentration, \f$ v \f$ is the flow velocity, 
		 \f$J_i\f$ is the molar flux, \f$D_{i,D}\f$ is the component-dust diffusion coefficient, 
         \f$D_{ij}\f$ is the pair diffusion coefficient. 

		 The second equation can be rewritten in the matrix form:
		 \f[ -\vec \nabla c_i = \sum_k A_{ik} \vec J_k,\;\;\; 
             A_{ik} \equiv \delta_{ik} \left(\frac{1}{D_{i,D}} + 
											 \Lambda_i + 
                                             \frac{c_i}{c_tD_{ik}} \right) - 
				          \frac{c_i}{c_tD_{ik}},\;\;\;
             \Lambda_i \equiv \sum_{j, j\neq i} \frac{c_j}{c_tD_{ij}}\f]
		 
		class parameters are related to the quation ones as follows 	 
		 \param cData corresponds to \f$c_i\f$
		 \param diffusionCoefficients corresponds to \f$D_{ij}\f$
		 \param velocity corresponds to \f$\vec v\f$	 
	*/
	class FDStefanMaxwell: public SingleKernelNM
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPAbstractDataWithGhostNodes Field;
			typedef acl::VectorOfElements Param;
			
		private:	
			std::vector<Data> cData;
			std::vector<Data> cInternalData;			

			Field efPhi;
			std::vector<Param> efCharge;

			Field velocity;

			const VectorTemplate* vectorTemplate;

			std::vector<std::vector<Param>> diffusionCoefficients;
			std::vector<Param> dustDiffusionCoefficients;

			virtual void init0();
			virtual void postProcessing();
		public:			
			FDStefanMaxwell();
			FDStefanMaxwell(Data c1,
			                Data c2,
			                const acl::VectorOfElements & dC, 
			                const VectorTemplate* vT);
			void setDiffusionCoefficient(acl::VectorOfElements d, 
			                             unsigned int i = 0,unsigned int j = 1);
			inline const Param & getDiffusionCoefficient(unsigned int i=0, unsigned int j=1) const;
			inline const Param & getDustDiffusionCoefficient(unsigned int i=0) const;
			void setDustDiffusionCoefficient(unsigned int i, const Param & dd);
			void setVectorTemplate(VectorTemplate* vT);
			inline const VectorTemplate* getVectorTemplate() const;
			void setElectricField(Field phi);
			Field getElectricField() const;
			inline const Param & getCharge(unsigned int i) const;
			void setCharge(unsigned int i, const Param & q);
			void setVelocity(Field v);			
			inline Field getVelocity();
			inline std::vector<Data> & getData();
			void addComponent(Data c, const Param & dC);
			void addComponent(Data c, const Param & dC, const Param & q);			
	};

	typedef std::shared_ptr<FDStefanMaxwell> SPFDStefanMaxwell;

	/**
		 \ingroup TransportProcesses
		 \ingroup NumMethods
		 \f[ \partial_t c_i= - \vec \nabla \cdot \vec J - \vec \nabla \cdot (\vec v c_i) \f]
		 \f[ -\nabla c_i = \sum_{j, i\neq j} \frac{c_j\vec J_i-c_i\vec J_j}{c_tD_{ij}} + \frac{\vec J_i}{D_{i,D}}\f]						
		 where \f$c_i\f$ is a molar concentration, \f$ v \f$ is the flow velocity, 
		 \f$J_i\f$ is the molar flux, \f$D_{i,D}\f$ is the component-dust diffusion coefficient, 
         \f$D_{ij}\f$ is the pair diffusion coefficient. 
		 
		parameters are related to the quation ones as follows 	 
		 \param c1 \f$c_1\f$
		 \param c2 \f$c_2\f$
		 \param diffusionCoeff corresponds to \f$D_{12}\f$
		 \param v velocity field
		 \param vt used VectorTemplate

	 */
	SPFDStefanMaxwell generateFDStefanMaxwell(SPDataWithGhostNodesACLData c1, 
	                                          SPDataWithGhostNodesACLData c2,
	                                          double diffustionCoeff,
	                                          SPAbstractDataWithGhostNodes v, 
	                                          const VectorTemplate* vt);
	
	/**
		 \ingroup TransportProcesses
		 \ingroup NumMethods
		 
		 \f[ \partial_t c_i= - \vec \nabla \cdot \vec J\f]
		 \f[ -\nabla c_i = \sum_{j, i\neq j} \frac{c_j\vec J_i-c_i\vec J_j}{c_tD_{ij}} + \frac{\vec J_i}{D_{i,D}}\f]						
		 where \f$c_i\f$ is a molar concentration, \f$ v \f$ is the flow velocity, 
		 \f$J_i\f$ is the molar flux, \f$D_{i,D}\f$ is the component-dust diffusion coefficient, 
         \f$D_{ij}\f$ is the pair diffusion coefficient. 
		 
		 \param c1 \f$c_1\f$
		 \param c2 \f$c_2\f$
		 \param diffusionCoeff corresponds to \f$D_{12}\f$
		 \param vt used VectorTemplate
	
	*/
	SPFDStefanMaxwell generateFDStefanMaxwell(SPDataWithGhostNodesACLData c1, 
	                                          SPDataWithGhostNodesACLData c2,
	                                          double diffustionCoeff,
	                                          const VectorTemplate* vt);

	class FDStefanMaxwellElectricField: public SingleKernelNM
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPAbstractDataWithGhostNodes Field;
			
		private:
			SPFDStefanMaxwell smSolver;
			Data phi;
			Data phiInternalData;
			Field phiS;
			const double stepFactor=1e-3;

			virtual void init0();
			virtual void postProcessing();
		public:			
			FDStefanMaxwellElectricField(SPFDStefanMaxwell sm, Data phi);
			void setPhiS(Field pS);
	};

	typedef std::shared_ptr<FDStefanMaxwellElectricField> SPFDStefanMaxwellElectricField;
	
// ------------------------- Implementation ------------------------

	inline FDStefanMaxwell::Field FDStefanMaxwell::getVelocity()
	{
		return velocity;
	}
	
	inline std::vector<FDStefanMaxwell::Data> & FDStefanMaxwell::getData()
	{
		return cData;
	}
	
	inline const VectorTemplate* FDStefanMaxwell::getVectorTemplate() const
	{
		return vectorTemplate;
	}

	inline const acl::VectorOfElements & 
		FDStefanMaxwell::getDiffusionCoefficient(unsigned int i, unsigned int j) const
	{
		return diffusionCoefficients[i][j];
	}

	inline const acl::VectorOfElements & 
		FDStefanMaxwell::getDustDiffusionCoefficient(unsigned int i) const
	{
		return dustDiffusionCoefficients[i];
	}

	inline const acl::VectorOfElements & 
		FDStefanMaxwell::getCharge(unsigned int i) const
	{
		return efCharge[i];
	}

	
} // asl
#endif // ASLFDADVECTIONDIFFUSION_H
