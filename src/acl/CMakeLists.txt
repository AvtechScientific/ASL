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