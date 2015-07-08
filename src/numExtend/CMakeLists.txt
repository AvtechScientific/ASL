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
INSTALL_SUBLIB(aslnumext aslnumext_PUBLIC_HEADERS)
