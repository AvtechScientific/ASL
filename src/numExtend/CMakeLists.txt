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
set_target_properties(aslExt PROPERTIES VERSION ${ASL_VERSION} SOVERSION ${ASL_VERSION_MAJOR})
INSTALL_SUBLIB(aslExt aslExt_PUBLIC_HEADERS)
