include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/acl)


# aslacloperators

set(aslacloperators_PUBLIC_HEADERS
	aclElementAssignmentSafe.h
	aclElementIfElse.h
	aclElementSum.h
	aclElementSubtraction.h
	aclElementSqrt.h
	aclElementSin.h
	aclElementProduct.h
	aclElementDivision.h
	aclElementGenericBinary.h
	aclGenericAtomicFunction.h
	aclElementGenericUnary.h
	aclElementCos.h
	aclElementSelect.h
	aclElementMad.h
	aclElementExcerpt.h
	aclElementParser.h
	aclElementSyncCopy.h
	aclOperatorUnary.h
	aclOperatorBinary.h
	aclOperatorTernary.h
	aclOperatorGeneric.h
	aclElementFor.h
	aclElementConvert.h
)

set(aslacloperators_SOURCES
	${aslacloperators_PUBLIC_HEADERS}
	aclElementAssignmentSafe.cxx
	aclElementIfElse.cxx
	aclElementSum.cxx
	aclElementSubtraction.cxx
	aclElementSqrt.cxx
	aclElementSin.cxx
	aclElementProduct.cxx
	aclElementDivision.cxx
	aclElementGenericBinary.cxx
	aclGenericAtomicFunction.cxx
	aclElementGenericUnary.cxx
	aclElementCos.cxx
	aclElementSelect.cxx
	aclElementMad.cxx
	aclElementExcerpt.cxx
	aclElementParser.cxx
	aclElementSyncCopy.cxx
	aclOperatorUnary.cxx
	aclOperatorBinary.cxx
	aclOperatorTernary.cxx
	aclOperatorGeneric.cxx
	aclElementFor.cxx
	aclElementConvert.cxx
)

add_library(aslacloperators ${aslacloperators_SOURCES})
target_link_libraries(aslacloperators PRIVATE ${OpenCL_LIBRARIES})
INSTALL_SUBLIB(aslacloperators aslacloperators_PUBLIC_HEADERS)
