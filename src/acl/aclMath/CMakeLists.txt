include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/acl)


# aslaclMath

set(aslaclMath_PUBLIC_HEADERS
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

set(aslaclMath_SOURCES
	${aslaclMath_PUBLIC_HEADERS}
	aclVectorOfElementsDef.cxx
	aclVectorOfElementsOperations.cxx
	aclMatrixOfElements.cxx
	aclComplexNumOfElements.cxx
	aclQuaternionOfElements.cxx
	aclReductionAlgGenerator.cxx
	aclBarycentric.cxx
	aclMathAlg.cxx
)

add_library(aslaclMath ${aslaclMath_SOURCES})
INSTALL_SUBLIB(aslaclMath aslaclMath_PUBLIC_HEADERS)
