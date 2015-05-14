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


#include "aclVariableSP.h"
#include "../aclUtilities.h"
#include "../aclTypesList.h"
#include "../aclHardware.h"


namespace acl
{


	template <> unsigned int VariableSP<cl_int>::id(0);
	template <> const string VariableSP<cl_int>::prefix("vsp_i");

	template <> unsigned int VariableSP<cl_uint>::id(0);
	template <> const string VariableSP<cl_uint>::prefix("vsp_ui");

	template <> unsigned int VariableSP<cl_float>::id(0);
	template <> const string VariableSP<cl_float>::prefix("vsp_f");

	template <> unsigned int VariableSP<cl_double>::id(0);
	template <> const string VariableSP<cl_double>::prefix("vsp_d");

	template <> unsigned int VariableSP<cl_long>::id(0);
	template <> const string VariableSP<cl_long>::prefix("vsp_l");

	template <typename T> string VariableSP<T>::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return TYPE[typeID] + " " + name;
	}
	#define BOOST_TT_rep_expression(r, data, T) \
	template string VariableSP<T>::getTypeSignature(const KernelConfiguration & kernelConfig) const;
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression
	
} // namespace acl
