#
# Advanced Simulation Library <http://asl.org.il>
# 
# Copyright 2015 Avtech Scientific <http://avtechscientific.com>
#
#
# This file is part of Advanced Simulation Library (ASL).
#
# ASL is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, version 3 of the License.
#
# ASL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with ASL. If not, see <http://www.gnu.org/licenses/>.
#


include(Operators/acloperators.cmake)
include(aclMath/aclmath.cmake)
include(DataTypes/acldatatypes.cmake)
include(Kernels/aclkernels.cmake)

include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/acl)


# aslacl

set(aslacl_PUBLIC_HEADERS
	${acloperators_PUBLIC_HEADERS}
	${aclmath_PUBLIC_HEADERS}
	${acldatatypes_PUBLIC_HEADERS}
	${aclkernels_PUBLIC_HEADERS}
	aclStdIncludes.h 
	aclTypesList.h 
	aclTypes.h 
	acl.h
	cl.hpp
	aclHardware.h
	aclUtilities.h
	aclElementBase.h
	aclGenerators.h
)

set(aslacl_SOURCES
	${acloperators_SOURCES}
	${aclmath_SOURCES}
	${acldatatypes_SOURCES}
	${aclkernels_SOURCES}
	acl.cxx 
	aclHardware.cxx
	aclUtilities.cxx 
	aclElementBase.cxx
	aclGenerators.cxx
)

add_library(aslacl ${aslacl_PUBLIC_HEADERS} ${aslacl_SOURCES})
target_link_libraries(aslacl PUBLIC aslcommon ${OpenCL_LIBRARIES})
INSTALL_SUBLIB(aslacl aslacl_PUBLIC_HEADERS)
