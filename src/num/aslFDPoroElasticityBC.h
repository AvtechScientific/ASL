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


#ifndef ASLFDPOROELASTICITYBC_H
#define ASLFDPOROELASTICITYBC_H

#include "aslBCond.h"
#include "acl/aclMath/aclVectorOfElementsDef.h"

namespace acl{
	class Kernel;
}

namespace asl
{
	class FDPoroElasticity;
	typedef std::shared_ptr<FDPoroElasticity> SPFDPoroElasticity;
	class DistanceFunction;
	typedef std::shared_ptr<DistanceFunction> SPDistanceFunction;
	class PositionFunction;
	typedef std::shared_ptr<PositionFunction> SPPositionFunction;
	
	/// Bondary condition corresponding to a rigid wall (\f$\vec u=0\f$ and \f$\nabla p=0\f$)
	/**	 
		 \ingroup ElasticityBC
	*/
	class BCRigidWallPoroElasticity:public BCond
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;			
			SPFDPoroElasticity num;
			acl::VectorOfElements value;
		public:
			BCRigidWallPoroElasticity(SPFDPoroElasticity nm);
			BCRigidWallPoroElasticity(SPFDPoroElasticity nm, acl::VectorOfElements v);
			~BCRigidWallPoroElasticity();
			virtual void execute();
			virtual void init();			
	};

	/// Bondary condition corresponding to a rigid wall (\f$\vec u=0\f$ and \f$\nabla p=0\f$)
	/**	 
		 \ingroup ElasticityBC
	*/
	class BCRigidWallDF: public BCondWithMap
	{		
		protected:
			std::unique_ptr<acl::Kernel> kernel;			
			SPFDPoroElasticity num;
			SPDistanceFunction rWall;
		public:
			BCRigidWallDF(SPFDPoroElasticity nm, 
			              SPDistanceFunction rw, 
			              SPAbstractDataWithGhostNodes map);
			~BCRigidWallDF();
			virtual void execute();
			virtual void init();			
	};
	
	void addBCRigidWall(std::vector<SPNumMethod> & bcList,
	                     SPFDPoroElasticity nm, 
	                     const std::vector<SlicesNames> & sl);
		
	/// Bondary condition set given values to pressure
	/**	 
		 \ingroup ElasticityBC
	*/
	void addBCZeroStress(std::vector<SPNumMethod> & bcList,
	                     SPFDPoroElasticity nm, 
	                     SPAbstractDataWithGhostNodes map);

	/// Bondary condition set given values to pressure
	/**	 
		 \ingroup ElasticityBC
		 \param bcList list of boundary conditions
		 \param nm corresponding numerical method
		 \param p pressure field
		 \param map computation map
	*/
	void addBCZeroStress(std::vector<SPNumMethod> & bcList,
	                     SPFDPoroElasticity nm,
	                     SPPositionFunction p,
	                     SPAbstractDataWithGhostNodes map);
	

	/**	 
		 \ingroup ElasticityBC
	*/
	void addBCRigidWallDF(std::vector<SPNumMethod> & bcList,
	                      SPFDPoroElasticity nm,
	                      SPDistanceFunction rw, 
			              SPAbstractDataWithGhostNodes map);

	/**	 
		 \ingroup ElasticityBC
	   considers \p rw as linear interpolation
	*/
	void addBCRigidWallDF(std::vector<SPNumMethod> & bcList,
	                      SPFDPoroElasticity nm,
	                      SPAbstractDataWithGhostNodes rw, 
			              SPAbstractDataWithGhostNodes map);
	
		 
} //asl

#endif //ASLFDPOROELASTICITYBC_H
