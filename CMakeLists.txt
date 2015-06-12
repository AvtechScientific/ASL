project(asl)

# After switching to 3.1 - remove FindOpenCL.cmake from the distribution
# since it is supported by the cmake itself.
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)

# Set version number
set(ASL_VERSION_MAJOR 0)
set(ASL_VERSION_MINOR 1)
set(ASL_VERSION_PATCH 3)
set(ASL_VERSION ${ASL_VERSION_MAJOR}.${ASL_VERSION_MINOR}.${ASL_VERSION_PATCH})

if (CMAKE_VERSION VERSION_LESS 3.1.0)
	set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
else()
	set(CMAKE_CXX_STANDARD 11)
endif()

# Add definitions
include(DefineOptions.cmake)
include(DefineFunctions.cmake)

# Configure Doxyfile
configure_file(
	"${PROJECT_SOURCE_DIR}/doc/Developer-Guide/Doxyfile.in"
	"${PROJECT_BINARY_DIR}/doc/Developer-Guide/Doxyfile"
)

# Enable supplied cmake Modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Find dependencies
find_package(OpenCL 1.1 REQUIRED)
# ToDo: check - system component might be not needed any longer on new versions of Boost
find_package(Boost 1.55 REQUIRED COMPONENTS program_options filesystem system)
find_package(VTK 6.1 COMPONENTS vtkRenderingCore vtkImagingCore vtkFiltersCore vtkIOCore vtkIOLegacy vtkIOXML vtkIOMINC vtkCommonCore vtkViewsCore vtkftgl vtksys vtkDICOMParser vtkexpat vtkzlib NO_MODULE)
include(${VTK_USE_FILE})

include_directories(SYSTEM ${OpenCL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

add_subdirectory(src)
add_subdirectory(test)

# Packaging
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION_MAJOR ${ASL_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${ASL_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${ASL_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")
include(CPack)
