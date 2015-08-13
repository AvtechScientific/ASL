include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/acl)

# aslnumext

set(aslnumext_PUBLIC_HEADERS
	aslFDAdvectionDiffusionExtended.h
)

set(aslnumext_SOURCES
	${aslnumext_PUBLIC_HEADERS}
	aslFDAdvectionDiffusionExtended.cxx
)

add_library(aslnumext ${aslnumext_SOURCES})
target_link_libraries(aslnumext PRIVATE asl aslacl aslcommon ${OpenCL_LIBRARIES})
INSTALL_SUBLIB(aslnumext aslnumext_PUBLIC_HEADERS)
