include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/acl)


# aslExt

set(aslExt_PUBLIC_HEADERS
	aslFDAdvectionDiffusionExtended.h
)

set(aslExt_SOURCES
	${aslExt_PUBLIC_HEADERS}
	aslFDAdvectionDiffusionExtended.cxx
)

add_library(aslExt ${aslExt_SOURCES})
INSTALL_SUBLIB(aslExt aslExt_PUBLIC_HEADERS)
