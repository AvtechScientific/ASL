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


#ifndef ASLFDADVECTIONDIFFUSIONBC_H
#define ASLFDADVECTIONDIFFUSIONBC_H

#include "aslBCond.h"

namespace acl
{
	class Kernel;
	typedef std::shared_ptr<Kernel> SPKernel;

}

namespace asl
{
	class FDAdvectionDiffusion;
	typedef std::shared_ptr<FDAdvectionDiffusion> SPFDAdvectionDiffusion;
	class AbstractDataWithGhostNodes;
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;

	
	/// Bondary condition that makes constant flux for pure diffusion \ingroup BoundaryConditions		 
	class BCConstantFluxMap:public BCondWithMap
	{		
		public:
			typedef SPAbstractDataWithGhostNodes Data;
			acl::SPKernel kernel;			
		protected:
			Data  data;			
			acl::VectorOfElements value; 
		public:
			BCConstantFluxMap(Data d,
			                  const acl::VectorOfElements & val, 
			                  Data map,
			                  const VectorTemplate *const t);
			~BCConstantFluxMap();
			virtual void execute();
			virtual void init();			
			void setValue(const acl::VectorOfElements & v);	
	};

	
	SPNumMethod generateBCConstantFlux(SPFDAdvectionDiffusion nm, 
	                                   double flux,
	                                   SPAbstractDataWithGhostNodes map);
		
} // asl
#endif // ASLFDADVECTIONDIFFUSIONBC_H
