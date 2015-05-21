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


#include "aclPrivateVariable.h"
#include "../aclUtilities.h"

namespace acl
{

	template <> unsigned int PrivateVariable<cl_int>::id(0);
	template <> const string PrivateVariable<cl_int>::prefix("pv_i");

	template <> unsigned int PrivateVariable<cl_uint>::id(0);
	template <> const string PrivateVariable<cl_uint>::prefix("pv_ui");

	template <> unsigned int PrivateVariable<cl_float>::id(0);
	template <> const string PrivateVariable<cl_float>::prefix("pv_f");

	template <> unsigned int PrivateVariable<cl_double>::id(0);
	template <> const string PrivateVariable<cl_double>::prefix("pv_d");

	template <> unsigned int PrivateVariable<cl_long>::id(0);
	template <> const string PrivateVariable<cl_long>::prefix("pv_l");
	
	
} // namespace acl
