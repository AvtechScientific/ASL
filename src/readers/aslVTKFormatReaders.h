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


#ifndef ASLVTKFORMATREADERS_H
#define ASLVTKFORMATREADERS_H


#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <acl/aclHardware.h>
#include <data/aslDataWrapper.h>
#include <aslGenerators.h>

using  namespace std;

namespace agl
{
	class TrianglesList;
	typedef std::shared_ptr<TrianglesList> SPTrianglesList;
}

namespace asl
{
	using namespace asl;

//    std::shared_ptr<DataWithGhostNodesACLData> SPDataWithGhostNodesACLData;

	/// Reads \p arrayNum of data from a file (detecting its format through its extension)
	/// to asl data and puts it in \p queue;
	/// returns asl data;
	/// Supported formats: .vtk .vti .mnc .dcm
	/// \ingroup IO
	SPDataWithGhostNodesACLData read(const string & fileName,
	                                 unsigned int arrayNum,
	                                 acl::CommandQueue queue = acl::hardware.defaultQueue);

	/// Reads \p arrayNum of data from a file (detecting its format through its extension)
	/// to asl data and puts it in \p queue;
	/// returns asl data;
	/// Supported formats: .vtp .stl
	/// \ingroup IO
	SPDataWithGhostNodesACLData readSurface(const string & fileName,
	                                        double dx,
	                                        acl::CommandQueue queue = acl::hardware.defaultQueue);

	SPDataWithGhostNodesACLData readSurface(const string & fileName,
	                                        Block & b,
	                                        acl::CommandQueue queue = acl::hardware.defaultQueue);

	// read surface with offsets retative to the surface outframe
	SPDataWithGhostNodesACLData readSurface(const string & fileName,
	                                        double dx,
	                                        double offset_X0, double offset_XE,
	                                        double offset_Y0, double offset_YE,
	                                        double offset_Z0, double offset_ZE,
	                                        acl::CommandQueue queue = acl::hardware.defaultQueue);
	
} // asl

#endif // ASLVTKFORMATREADERS_H

