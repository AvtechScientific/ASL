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


#include "aslGenerators.h"
#include "acl/aclGenerators.h"
#include "acl/aclMath/aclMathAlg.h"
#include <data/aslDataWithGhostNodes.h>
#include "acl/aclHardware.h"

#include "acl/aclTypesList.h"


namespace asl
{

	SPDataWrapperACL generateDataContainer_SP(const Block &b, const acl::VectorOfElements & a)
	{
		auto aa(std::make_shared<DataWrapperACL>(b));
		aa->setContainer (a);
		return aa;
	}
	
	SPDataWithGhostNodesACL generateDataContainer_SP(const Block &b, const acl::VectorOfElements & a, unsigned int gN)
	{
		auto aa(std::make_shared<DataWithGhostNodesACL>(b,gN));
		aa->setContainer (a);
		return aa;
	}

	SPDataWithGhostNodesACLData generateDataContainerACL_SP(const Block &b, 
	                                                        acl::TypeID t, 
	                                                        unsigned int n, 
	                                                        unsigned int gN,
	                                                        acl::CommandQueue queue)
	{
		auto a(std::make_shared<DataWithGhostNodesACLData>(b,gN));
		a->setContainer(acl::generateVEData(productOfElements(a->getBlock().getSize()),
		                                    t,
		                                    n,
		                                    queue
		                                    ));
		return a;
	}

	SPDataWithGhostNodesACLData generateDataContainerACL_SP(const Block &b, 
	                                                        acl::TypeID t, 
	                                                        unsigned int n, 
	                                                        unsigned int gN)
	{
		return generateDataContainerACL_SP(b, t, n, gN, acl::hardware.defaultQueue);
	}
	

	SPAbstractDataWithGhostNodes generateDCFullSafe(SPAbstractDataWithGhostNodes d)
	{
		auto & b(d->getBlock());

		auto a(std::make_shared<DataWithGhostNodesACL>(b,0));
		a->setContainer(acl::generateVEOutOfBoundarySafe(d->getEContainer()));
		
		return a;
	}

	SPAbstractDataWithGhostNodes generateDCFullSafe(SPAbstractDataWithGhostNodes d, double outVal)
	{
		auto & b(d->getBlock());

		auto a(std::make_shared<DataWithGhostNodesACL>(b,0));
		a->setContainer(acl::generateVEOutOfBoundarySafe(d->getEContainer(), 
		                                                 acl::generateVEConstant(outVal)));
		
		return a;
	}

	template <typename T> SPDataWrapperACLData 
		generateDataContainerACL_SP(const Block &b, unsigned int n)
	{
		auto a(std::make_shared<DataWrapperACLData>(b));
		a->setContainer(acl::generateVEData<T>(productOfElements(b.getSize()),n));
		return a;
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template SPDataWrapperACLData generateDataContainerACL_SP<t>(const Block &b, unsigned int n);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
	
	template <typename T> SPDataWithGhostNodesACLData 
		generateDataContainerACL_SP(const Block &b,
		                            unsigned int n, 
	                                unsigned int gN,
	                                acl::CommandQueue queue)
	{
		auto a(std::make_shared<DataWithGhostNodesACLData>(b,gN));
		a->setContainer(acl::VectorOfElementsData(n,
		                                          productOfElements(a->getBlock().getSize()),
		                                          T(0), 
		                                          queue));
		return a;
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template SPDataWithGhostNodesACLData generateDataContainerACL_SP<t>(const Block &b, unsigned int n, unsigned int gN, acl::CommandQueue queue);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression

	template <typename T> SPDataWithGhostNodesACLData 
		generateDataContainerACL_SP(const Block &b,
		                            unsigned int n, 
	                                unsigned int gN)
	{
		return generateDataContainerACL_SP<T>(b, n, gN, acl::hardware.defaultQueue);
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template SPDataWithGhostNodesACLData generateDataContainerACL_SP<t>(const Block &b, unsigned int n, unsigned int gN);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
		
	template <typename T> SPDataWithGhostNodesACL 
		generateDataContainerConst_SP(const Block &b,T a, unsigned int gN)
	{
		auto aa(std::make_shared<DataWithGhostNodesACL>(b,gN));
		aa->setContainer(acl::generateVEConstant(a));
		return aa;
	}

	#define BOOST_TT_rep_expression(r, data, t) \
		template SPDataWithGhostNodesACL generateDataContainerConst_SP(const Block &b,t a, unsigned int gN);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
	#define BOOST_TT_rep_expression(r, data, t) \
		template SPDataWithGhostNodesACL generateDataContainerConst_SP(const Block &b, AVec<t> a, unsigned int gN);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
	
	
} // namespace asl
