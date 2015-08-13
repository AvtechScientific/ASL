# aslmath

set(aslmath_PUBLIC_HEADERS
	math/aslVectors.h
	math/aslVectorsDynamicLength.h
	math/aslVectorsDynamicLengthOperations.h
	math/aslMatrices.h
	math/aslBarycentric.h
	math/aslInterpolation.h
	math/aslTemplates.h
	math/aslTemplatesExtras.h
	math/aslProbeTemplates.h
	math/aslTemplateVE.h
	math/aslTemplateVEExtras.h
	math/aslIndex2Position.h
	math/aslDistanceFunction.h
	math/aslDistanceFunctionAlg.h
	math/aslPositionFunction.h
)

set(aslmath_SOURCES
	math/aslMatrices.cxx
	math/aslBarycentric.cxx
	math/aslTemplates.cxx
	math/aslTemplatesExtras.cxx
	math/aslProbeTemplates.cxx
	math/aslTemplateVE.cxx
	math/aslTemplateVEExtras.cxx
	math/aslIndex2Position.cxx
	math/aslDistanceFunction.cxx
	math/aslDistanceFunctionAlg.cxx
	math/aslPositionFunction.cxx
)