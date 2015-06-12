add_subdirectory(Operators)
add_subdirectory(aclMath)

include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/acl)


# aslacl

set(aslacl_PUBLIC_HEADERS
	aclStdIncludes.h 
	aclTypesList.h 
	aclTypes.h 
	acl.h 
	aclHardware.h
	aclUtilities.h
	DataTypes/aclConstant.h
	DataTypes/aclIndex.h 
	DataTypes/aclGroupID.h
	DataTypes/aclIndexExt.h 
	DataTypes/aclPrivateVariable.h 
	DataTypes/aclPrivateArray.h
	DataTypes/aclLocalArray.h
	DataTypes/aclVariable.h 
	DataTypes/aclVariableReference.h 
	DataTypes/aclVariableSP.h
	DataTypes/aclMemBlock.h
	DataTypes/aclArray.h
	DataTypes/aclSubvector.h
	Kernels/aclExpressionContainer.h
	Kernels/aclKernel.h
	Kernels/aclKernelConfiguration.h
	Kernels/aclKernelConfigurationTemplates.h
	Kernels/aclKernelMerger.h
	aclElementBase.h
	aclGenerators.h
)

set(aslacl_SOURCES
	${aslacl_PUBLIC_HEADERS}
	acl.cxx 
	aclHardware.cxx
	aclUtilities.cxx 
	DataTypes/aclConstant.cxx
	DataTypes/aclIndex.cxx
	DataTypes/aclGroupID.cxx
	DataTypes/aclIndexExt.cxx
	DataTypes/aclPrivateVariable.cxx 
	DataTypes/aclPrivateArray.cxx
	DataTypes/aclLocalArray.cxx
	DataTypes/aclVariable.cxx 
	DataTypes/aclVariableReference.cxx 
	DataTypes/aclVariableSP.cxx 
	DataTypes/aclMemBlock.cxx
	DataTypes/aclArray.cxx
	DataTypes/aclSubvector.cxx
	Kernels/aclExpressionContainer.cxx
	Kernels/aclKernel.cxx
	Kernels/aclKernelConfiguration.cxx
	Kernels/aclKernelConfigurationTemplates.cxx
	Kernels/aclKernelMerger.cxx
	aclElementBase.cxx
	aclGenerators.cxx
)

add_library(aslacl ${aslacl_SOURCES})
set_target_properties(aslacl PROPERTIES VERSION ${ASL_VERSION} SOVERSION ${ASL_VERSION_MAJOR})
target_link_libraries(aslacl aslCommon aslaclOperators aslaclMath ${OpenCL_LIBRARIES})
INSTALL_SUBLIB(aslacl aslacl_PUBLIC_HEADERS)
