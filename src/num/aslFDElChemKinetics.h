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


#ifndef ASLFDELCHEMKINETICS_H
#define ASLFDELCHEMKINETICS_H

#include "aslSingleKernelNM.h"
#include "acl/aclMath/aclVectorOfElementsDef.h"


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
	
	/// Numerical method which computes electrode reactions
	/** 
		 !!!!!!!!!!!!need corrections

		 
		 \ingroup ChemicalReactions
		 \ingroup NumMethods

		 Description of an electrochemical kinetic reaction on an electrode:
		 \f[ \ce{A <-> A^{n+} + n e-} \f]
		 The reaction kinetics is described by the Butler-Volmer equation:
		 \f[ j = j_0 c_A \cdot \left( e^{\frac { \alpha_a nF \eta} {RT}} - e^{ - { \frac { \alpha_c nF \eta} {RT}} } \right), \f]
		 where
			 
		- \f$ j_0 \f$ is an exchange current density, \f$ A/m^2 \f$, note the algorithm assumes \f$ mol/m^2/s \f$  
		- \f$ T \f$ is absolute temperature, \f$K\f$ 
		- \f$ n \f$ is number of electrons involved in the electrode reaction
		- \f$ F \f$ is the Faraday constant, \f$ 9.64853399 \cdot 10^4 C mol^{-1} \f$
		- \f$ R \f$ is universal gas constant, \f$ 8.3144621 J K^{-1} mol^{-1} \f$  
		- \f$\alpha_c\f$ is so-called cathodic charge transfer coefficient
		- \f$\alpha_a\f$ is so-called anodic charge transfer coefficient, dimensionless
		- \f$\eta\f$ is activation overpotential defined as: \f$\eta = \phi - \phi_{eq}\f$		 
		- \f$E\f$ is electrode potential, \f$ V \f$
		- \f$E_{eq}\f$ equilibrium potential, \f$ V \f$
		- \f$c_A\f$ concentration of the component \f$\ce{A}\f$

		The equation describes chemical reaction namely it can be used together 
        with an appropriate transport equations and equation for the elecrostatic potential:
		\f[ c_A^R = c_A - j \f]
		\f[ c_{A^{n+}}^R = c_{A^{n+}} + j \f]
		\f[ \vec\nabla \kappa \vec\nabla\phi + S_\phi = 0 \f]

		class parameters are related to the quation ones as follows 	 
		 \param a0 molar concentration of the component before reaction
		 \param aI molar concentration of the ionized component
		 \param efSPhi source term in the equation for electrostatic potential
		 \param phi electrostatic potential, \f$ \phi \f$	
		 \param j0 corresponds to \f$ j_0 \f$, note it should have the same units as concentrations
		 \param phi0 corresponds to \f$ \phi_{eq} \f$
		 \param alphaA corresponds to \f$ \frac { \alpha_a F} {RT} \f$ 
		 \param alphaC corresponds to \f$ \frac { \alpha_c F} {RT} \f$
		 \param n corresponds to \f$ n \f$
			
	*/
	class FDBVKinetics: public SingleKernelMapNM
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPAbstractDataWithGhostNodes Field;
			typedef acl::VectorOfElements Param;
			
		private:
			acl::SPKernel kernelJ;
			std::vector<Data> aI;			
			Data efSPhi;
			Field phi;

			vector<double> nI;
			Param j0;
			Param beta;
			double n;

			virtual void init0();
		public:			
			FDBVKinetics(Data a0,
			             double n0,
			             Data aI,
			             double nI,
			             Data phiS,
			             Field phi,
			             const Param & j0,
			             const Param & beta,
			             double n);
			void setElectricFieldSource(Field phi);
			Field getElectricFieldSource() const;
			inline Data & getAI(unsigned int i=0);
			void addAI(Data ai, double ni);			
			void executeJ();
	};

	typedef std::shared_ptr<FDBVKinetics> SPFDBVKinetics;


	SPFDBVKinetics generateFDBVKinetics(SPDataWithGhostNodesACLData a0,
	                                    double n0,
	                                    SPDataWithGhostNodesACLData aI,
	                                    double nI,
	                                    SPDataWithGhostNodesACLData phiS,
	                                    SPAbstractDataWithGhostNodes phi,
	                                    double j0,
	                                    double beta,
	                                    double n);
	

	
// ------------------------- Implementation ------------------------
	
	inline FDBVKinetics::Data & FDBVKinetics::getAI(unsigned int i)
	{
		return aI[i];
	}
		
} // asl
#endif // ASLFDADVECTIONDIFFUSION_H
