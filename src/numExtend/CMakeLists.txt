include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/acl)

# aslnumext

set(aslnumext_PUBLIC_HEADERS
	aslFDAdvectionDiffusionExtended.h
)

set(aslnumext_SOURCES
	aslFDAdvectionDiffusionExtended.cxx
)

add_library(aslnumext ${aslnumext_PUBLIC_HEADERS} ${aslnumext_SOURCES})
target_link_libraries(aslnumext PRIVATE asl aslacl aslcommon ${OpenCL_LIBRARIES})
INSTALL_SUBLIB(aslnumext aslnumext_PUBLIC_HEADERS)
