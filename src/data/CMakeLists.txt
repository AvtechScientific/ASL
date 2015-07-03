include_directories(${CMAKE_SOURCE_DIR}/src)


# asldata

set(asldata_PUBLIC_HEADERS
	aslMemElement.h
	aslDataUtilities.h
	aslBlocks.h
	aslDataWrapper.h
	aslDataWithGhostNodes.h
	aslProbe.h
)

set(asldata_SOURCES
	${asldata_PUBLIC_HEADERS}
	aslMemElement.cxx
	aslDataUtilities.cxx
	aslBlocks.cxx
	aslDataWrapper.cxx
	aslDataWithGhostNodes.cxx
	aslProbe.cxx
)

add_library(asldata ${asldata_SOURCES})
INSTALL_SUBLIB(asldata asldata_PUBLIC_HEADERS)
