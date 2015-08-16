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


#ifndef ACLTYPES_H
#define ACLTYPES_H

namespace acl
{
	enum Extension
	{
		CL_KHR_FP64,
		CL_KHR_INT64_BASE_ATOMICS,
		CL_KHR_INT64_EXTENDED_ATOMICS,
		CL_KHR_GL_SHARING,
	};


	enum TypeID
	{
		TYPE_INT,
		TYPE_UINT,
		TYPE_FLOAT,
		TYPE_DOUBLE,
		TYPE_LONG,
	};
}

#endif // ACLTYPES_H