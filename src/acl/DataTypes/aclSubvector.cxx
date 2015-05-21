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


#include "aclSubvector.h"
#include "../../aslUtilities.h"

namespace acl
{

	
	template <> unsigned int Subvector<cl_int>::id(0);
	template <> const string Subvector<cl_int>::prefix("svi");

	template <> unsigned int Subvector<cl_uint>::id(0);
	template <> const string Subvector<cl_uint>::prefix("svui");
	
	template <> unsigned int Subvector<cl_float>::id(0);
	template <> const string Subvector<cl_float>::prefix("svf");

	template <> unsigned int Subvector<cl_double>::id(0);
	template <> const string Subvector<cl_double>::prefix("svd");

	template <> unsigned int Subvector<cl_long>::id(0);
	template <> const string Subvector<cl_long>::prefix("svl");

}
