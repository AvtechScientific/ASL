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


#include "aclConstant.h"
#include "../aclUtilities.h"

#include <math.h>
#include "aclTypesList.h"


namespace acl
{

	template <typename T> Constant<T>::Constant(T c):
		ElementBase(false, 0, typeToTypeID<T>()),
		value(c)
	{
		stringstream s;
		if( value < 0 )
			s << "(" << value << ")";
		else
			s << value;
		valueStr = s.str();
	}

	template <> Constant<cl_double>::Constant(cl_double c):
		ElementBase(false, 0, typeToTypeID<cl_double>()),
		value(c)
	{
		stringstream s;
		double intp(0.);
		if (modf(value, &intp) == 0. && std::fabs(value) < 10000)
			s << value << ".";
		else
		{
			s.precision(20);
			s << value;
		}
		valueStr = s.str();
	}

	template <> Constant<cl_float>::Constant(cl_float c):
		ElementBase(false, 0, typeToTypeID<cl_float>()),
		value(c)
	{
		stringstream s;
		float intp(0.);
		if (modf(value, &intp) == 0. && std::fabs(value) < 10000)
			s << value << ".";
		else
		{
			s.precision(10);
			s << value;
		}
		valueStr = s.str();
	}

	#define BOOST_TT_rep_expression(r,data,t) \
		template Constant<t>::Constant(t c);
	BOOST_PP_SEQ_FOR_EACH(BOOST_TT_rep_expression, ~, BOOST_TT_acl_types)
	#undef BOOST_TT_rep_expression	
		
} // namespace acl
