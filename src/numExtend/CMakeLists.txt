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
target_link_libraries(aslnumext PRIVATE ${OpenCL_LIBRARIES} aslacl aslcommon asl)
INSTALL_SUBLIB(aslnumext aslnumext_PUBLIC_HEADERS)
