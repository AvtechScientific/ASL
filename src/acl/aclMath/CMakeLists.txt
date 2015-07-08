include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/acl)

# aslaclmath

set(aslaclmath_PUBLIC_HEADERS
	aclVectorOfElements.h
	aclVectorOfElementsDef.h
	aclVectorOfElementsOperations.h
	aclMatrixOfElements.h
	aclComplexNumOfElements.h
	aclQuaternionOfElements.h
	aclReductionAlgGenerator.h
	aclBarycentric.h
	aclMathAlg.h
)

set(aslaclmath_SOURCES
	${aslaclmath_PUBLIC_HEADERS}
	aclVectorOfElementsDef.cxx
	aclVectorOfElementsOperations.cxx
	aclMatrixOfElements.cxx
	aclComplexNumOfElements.cxx
	aclQuaternionOfElements.cxx
	aclReductionAlgGenerator.cxx
	aclBarycentric.cxx
	aclMathAlg.cxx
)

add_library(aslaclmath ${aslaclmath_SOURCES})
INSTALL_SUBLIB(aslaclmath aslaclmath_PUBLIC_HEADERS)
