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


#ifndef TEMPLATEVEEXTRAS_H_INCLUDED
#define TEMPLATEVEEXTRAS_H_INCLUDED

#include <acl/aclMath/aclVectorOfElements.h>
#include "data/aslDataWithGhostNodes.h"
#include "aslTemplates.h"

/**
 \defgroup DifferentialOpperators Differential Operrators
 */

namespace asl
{
	/// returns VectorOfElements containing values in cell corners \ingroup DifferentialOpperators 
	acl::VectorOfElements cellValues(const TemplateVE & a, unsigned int iEl);

	/// differential operator \f$ \vec\nabla a \f$ within a cell \ingroup DifferentialOpperators
	acl::VectorOfElements gradient(const TemplateVE & a, unsigned int iEl);

	/// differential operator \f$ \vec\nabla a \f$ within all cells \ingroup DifferentialOpperators
	/**
		The function generates code that computes the gradient values in each cell. The
		 gradient values are written in \p values which are private variables. The corresponding 
		 amount of variables are automaticaly generated.
	*/
	vector<acl::Element> gcGradientAllCells(const TemplateVE & a, 
	                                        vector<acl::VectorOfElements> & values);
	
		
}// asl

#endif // TEMPLATEVEEXTRAS_H_INCLUDED
