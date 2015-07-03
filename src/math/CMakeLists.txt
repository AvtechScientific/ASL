
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/acl)


# aslmath

set(aslmath_PUBLIC_HEADERS
	aslVectors.h
	aslVectorsDynamicLength.h
	aslVectorsDynamicLengthOperations.h
	aslMatrices.h
	aslBarycentric.h
	aslInterpolation.h
	aslTemplates.h
	aslTemplatesExtras.h
	aslProbeTemplates.h
	aslTemplateVE.h
	aslTemplateVEExtras.h
	aslIndex2Position.h
	aslDistanceFunction.h
	aslDistanceFunctionAlg.h
	aslPositionFunction.h
)

set(aslmath_SOURCES
	${aslmath_PUBLIC_HEADERS}
	aslMatrices.cxx
	aslBarycentric.cxx
	aslTemplates.cxx
	aslTemplatesExtras.cxx
	aslProbeTemplates.cxx
	aslTemplateVE.cxx
	aslTemplateVEExtras.cxx
	aslIndex2Position.cxx
	aslDistanceFunction.cxx
	aslDistanceFunctionAlg.cxx
	aslPositionFunction.cxx
)

add_library(aslmath ${aslmath_SOURCES})
INSTALL_SUBLIB(aslmath aslmath_PUBLIC_HEADERS)