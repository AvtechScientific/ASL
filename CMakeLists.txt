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


project(ASL)

# After switching to 3.1 - remove FindOpenCL.cmake from the distribution
# since it is supported by the cmake itself.
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)

# Set version number
set(ASL_VERSION_MAJOR 0)
set(ASL_VERSION_MINOR 1)
set(ASL_VERSION_PATCH 5)

set(ASL_VERSION ${ASL_VERSION_MAJOR}.${ASL_VERSION_MINOR}.${ASL_VERSION_PATCH})

if (CMAKE_VERSION VERSION_LESS 3.1.0)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -std=c++11")
else()
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

# Enable supplied cmake includes and modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/" "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Add ASL cmake definitions
include(ASLPaths)
include(ASLBuildOptions)
include(ASLBuildFunctions)

# Find dependencies
find_package(OpenCL 1.1 REQUIRED)
# ToDo: check - system component might be not needed any longer on new versions of Boost
find_package(Boost 1.55 REQUIRED COMPONENTS program_options filesystem system)
find_package(VTK 6.1 COMPONENTS vtkRenderingCore vtkImagingCore vtkFiltersCore vtkIOCore vtkIOLegacy vtkIOXML vtkIOMINC vtkCommonCore vtkViewsCore vtkftgl vtksys vtkDICOMParser vtkexpat vtkzlib NO_MODULE)
include(${VTK_USE_FILE})

include_directories(SYSTEM ${OpenCL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

add_subdirectory(src)
add_subdirectory(utilities)

if (WITH_EXAMPLES)
	add_subdirectory(examples)
endif()

if (WITH_API_DOC)
	add_subdirectory(doc)
endif()

if (WITH_TESTS)
	# Test labels:
	# DoublePrecision - tests that use double precision (will fail if device does not support it)
	# IO - Input/Output tests (may require input files)
	# Performance - performance tests
	#
	# Run them e.g.: `ctest [-C Release|Debug] -L Performance`
	# Exclude some of them: `ctest -LE DoublePrecision -LE Performance`
	enable_testing()
	add_subdirectory(test)
endif()

# Packaging
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION_MAJOR ${ASL_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${ASL_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${ASL_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Advanced Simulation Library is a free and open source OpenCL-based multiphysics simulation software package.")
set(CPACK_PACKAGE_DESCRIPTION "Advanced Simulation Library (ASL) is a free and open source multiphysics simulation software package. Its computational engine is based, among others, on the Lattice Boltzmann Methods and is written in OpenCL which enable extraordinarily efficient deployment on a variety of massively parallel architectures, ranging from inexpensive FPGAs, DSPs and GPUs up to heterogeneous clusters and supercomputers. The engine is hidden entirely behind C++ classes, so that no OpenCL knowledge is required from application programmers. ASL can be utilized to model various coupled physical and chemical phenomena and employed in a multitude of fields: computational fluid dynamics, virtual sensing, industrial process data validation and reconciliation, image-guided surgery, computer-aided engineering, high-performance scientific computing, etc..")
set(CPACK_PACKAGE_VENDOR "Avtech Scientific")
set(CPACK_PACKAGE_CONTACT "http://avtechscientific.com/contact")

## .deb-specific
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "http://asl.org.il")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_SUGGESTS "paraview")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}\n${CPACK_PACKAGE_DESCRIPTION}")

## .rpm-specific
set(CPACK_RPM_PACKAGE_URL "http://asl.org.il")

## NSIS-specific
set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\asl.org.il")

include(CPack)


# Config (after package definitions, because CPACK_PACKAGE_DESCRIPTION_SUMMARY
# is used for pkgconfig.

## ASLConfig.cmake
include(CMakePackageConfigHelpers)

# ASL conforms to Semantic Versioning scheme <http://semver.org/>
# Switch to COMPATIBILITY SameMajorVersion on ASL_VERSION 1.0.0
# Till then use COMPATIBILITY ExactVersion 
write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/ASLConfigVersion.cmake"
	VERSION ${ASL_VERSION}
	COMPATIBILITY ExactVersion
)

export(EXPORT ASLTargets
	FILE "${CMAKE_CURRENT_BINARY_DIR}/ASLTargets.cmake"
	NAMESPACE ASL::
)

set(ConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/ASL)

configure_package_config_file(
	${CMAKE_SOURCE_DIR}/cmake/ASLConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/ASLConfig.cmake
	INSTALL_DESTINATION ${ConfigPackageLocation}
)

install(EXPORT ASLTargets
	FILE
	ASLTargets.cmake
	NAMESPACE
	ASL::
	DESTINATION
	${ConfigPackageLocation}
)

install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/ASLConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/ASLConfigVersion.cmake"
	DESTINATION
	${ConfigPackageLocation}
	COMPONENT
	Devel
)


## ASL.pc

### Read the composed list of all sublibs from the ASL_SUBLIBS_GLOBAL_PROPERTY
get_property(ASL_SUBLIBS GLOBAL PROPERTY ASL_SUBLIBS_GLOBAL_PROPERTY)

configure_file(
	${CMAKE_SOURCE_DIR}/cmake/ASL.pc.in
	${CMAKE_CURRENT_BINARY_DIR}/ASL.pc
	@ONLY
)

install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/ASL.pc
	DESTINATION
	${CMAKE_INSTALL_LIBDIR}/pkgconfig
)