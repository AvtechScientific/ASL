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


#ifndef ACLREDUCTIONALGGENERATUR_H
#define ACLREDUCTIONALGGENERATUR_H

#include "aclVectorOfElementsDef.h"
#include <math/aslVectors.h>
#include <utilities/aslUValue.h>


namespace acl
{
	class Kernel;
	enum ReductionOperatorType {ROT_SUM, ROT_PRODUCT, ROT_MINIMUM, ROT_MAXIMUM};
		
	/// The class generates code corresponding to a reduction operation of elements 
	/**
		 \ingroup ComplexDataTypes
	 */
	template <typename ResType, enum ReductionOperatorType Op> class ReductionAlgGenerator
	{
		private:
			VectorOfElements ve;
			unsigned int nGroups;
			unsigned int nLocal;
		public:
		    asl::UValue<asl::AVec<ResType>> res;
		private:
			std::vector<std::vector<ResType>> groupRes;
			VectorOfElementsData groupResACL;

			shared_ptr<Kernel> kernel;
			
		public:	
			ReductionAlgGenerator(VectorOfElements v);
			void compute();	
			void generateAlg(Kernel & k);
			void generateAlg();
	};

	template <typename ResType> std::shared_ptr<ReductionAlgGenerator<ResType, ROT_SUM>> 
		generateSumAlg(VectorOfElements v);

	template <typename ResType> std::shared_ptr<ReductionAlgGenerator<ResType, ROT_MINIMUM>> 
		generateMinAlg(VectorOfElements v);

	template <typename ResType> std::shared_ptr<ReductionAlgGenerator<ResType, ROT_MAXIMUM>> 
		generateMaxAlg(VectorOfElements v);

	template <typename ResType> std::shared_ptr<ReductionAlgGenerator<ResType, ROT_PRODUCT>> 
		generateProductAlg(VectorOfElements v);
	

}  //namespace acl

#endif // ACLREDUCTIONALGGENERATUR_H
