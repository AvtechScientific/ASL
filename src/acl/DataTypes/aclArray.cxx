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


#include "aclArray.h"
#include "../../aslUtilities.h"

namespace acl
{
	template <> unsigned int Array<cl_int>::id(0);
	template <> const string Array<cl_int>::prefix("a_i");

	template <> unsigned int Array<cl_uint>::id(0);
	template <> const string Array<cl_uint>::prefix("a_ui");

	template <> unsigned int Array<cl_float>::id(0);
	template <> const string Array<cl_float>::prefix("a_f");

	template <> unsigned int Array<cl_double>::id(0);
	template <> const string Array<cl_double>::prefix("a_d");

	template <> unsigned int Array<cl_long>::id(0);
	template <> const string Array<cl_long>::prefix("a_l");
	
}
