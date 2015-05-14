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


#ifndef ASLSINGLEKERNELNM_H
#define ASLSINGLEKERNELNM_H

#include "aslNumMethod.h"

namespace acl
{
	class ExpressionContainer;
	class Kernel;
	class KernelConfiguration;
	typedef std::shared_ptr<Kernel> SPKernel;
}


namespace asl
{
	class AbstractDataWithGhostNodes;
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;

	
	/// Virtual class describes general interface for Numerical methods whith single Kernel;
	/// \ingroup Numerics
	class SingleKernelNM:public NumMethod
	{
		protected: 
			acl::SPKernel kernel;
			/// the function executed before kernel->compute()
			virtual void preProcessing();
			/// the function executed after kernel->compute()
			virtual void postProcessing();
			/// full initialisation but without kernel->setup()
			virtual void init0()=0;

			SingleKernelNM(const acl::KernelConfiguration & kernelCongig);
		public:
			virtual void execute();
			/// Builds the necesery internal data and kernels
			virtual void init();
			virtual ~SingleKernelNM();

			friend class NumMethodsMerger;
	};

	class SingleKernelMapNM:public SingleKernelNM
	{
		public:
			typedef SPAbstractDataWithGhostNodes Field;
		private:
			Field map;
		protected:
			SingleKernelMapNM(const acl::KernelConfiguration & kernelCongig);
			SingleKernelMapNM(Field map, const acl::KernelConfiguration & kernelCongig);
			void initMapInfrastructure(acl::ExpressionContainer & k);	
		public:
			void setMap(Field map);
			~SingleKernelMapNM();
	};
	


//---------------------------- Implementation ---------------------------


	
}	//asl

#endif //ASLSINGLEKERNELNM_H
