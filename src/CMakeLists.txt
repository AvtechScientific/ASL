# ToDo: recheck here and everywhere whether all listed
# _PUBLIC_HEADERS are indeed needed and if not move them
# to the _SOURCES (or _PRIVATE_HEADERS?)

add_subdirectory(acl)
add_subdirectory(data)
add_subdirectory(math)
add_subdirectory(num)
add_subdirectory(numExtend)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/acl)


# aslCommon

set(aslCommon_PUBLIC_HEADERS
	utilities/aslParametersManager.h
	utilities/aslTimeStamp.h
	aslUtilities.h
)

set(aslCommon_SOURCES
	${aslCommon_PUBLIC_HEADERS}
	utilities/aslParametersManager.cxx
	utilities/aslTimeStamp.cxx
	aslUtilities.cxx
)

add_library(aslCommon ${aslCommon_SOURCES})
set_target_properties(aslCommon PROPERTIES VERSION ${ASL_VERSION} SOVERSION ${ASL_VERSION_MAJOR})
target_link_libraries(aslCommon ${Boost_LIBRARIES})
INSTALL_SUBLIB(aslCommon aslCommon_PUBLIC_HEADERS)


# asl

set(asl_PUBLIC_HEADERS
	aslGenerators.h
	utilities/aslTimer.h
	utilities/aslUValue.h
	utilities/aslGlobalSpace.h
	utilities/aslSmartPtrUtils.h
	writers/aslABDFormat.h
	writers/aslWriter.h
	aslDataInc.h
	aslGeomInc.h
)

set(asl_SOURCES
	${asl_PUBLIC_HEADERS}
	aslGenerators.cxx
	writers/aslABDFormat.cxx
	writers/aslWriter.cxx
	aslDataInc.cxx
	aslGeomInc.cxx
)

add_library(asl ${asl_SOURCES})
set_target_properties(asl PROPERTIES VERSION ${ASL_VERSION} SOVERSION ${ASL_VERSION_MAJOR})
target_link_libraries(asl aslCommon aslacl asldata aslmath ${Boost_LIBRARIES})
INSTALL_SUBLIB(asl asl_PUBLIC_HEADERS)


# aslvtk

set(aslvtk_PUBLIC_HEADERS
	utilities/aslVTKCasters.h
	utilities/aslVTKDataGenerators.h
	writers/aslVTKFormatWriters.h
	readers/aslVTKFormatReaders.h
)

set(aslvtk_SOURCES
	${aslvtk_PUBLIC_HEADERS}
	utilities/aslVTKCasters.cxx
	utilities/aslVTKDataGenerators.cxx
	writers/aslVTKFormatWriters.cxx
	readers/aslVTKFormatReaders.cxx
)

add_library(aslvtk ${aslvtk_SOURCES})
set_target_properties(aslvtk PROPERTIES VERSION ${ASL_VERSION} SOVERSION ${ASL_VERSION_MAJOR})
target_link_libraries(aslvtk ${VTK_LIBRARIES} ${OpenCL_LIBRARIES})
INSTALL_SUBLIB(aslvtk aslvtk_PUBLIC_HEADERS)


# aslmatio

if (WITH_MATIO)
	find_package(MATIO 1.5.2 REQUIRED)
	include_directories(SYSTEM ${MATIO_INCLUDE_DIRS} ${OpenCL_LIBRARIES})
	set(CMAKE_CXX_FLAGS "-fpermissive ${CMAKE_CXX_FLAGS}")

	set(aslmatio_PUBLIC_HEADERS
		utilities/aslMATLABCasters.h
		writers/aslMATFormat.h
	)

	set(aslmatio_SOURCES
		${aslmatio_PUBLIC_HEADERS}
		utilities/aslMATLABCasters.cxx
		writers/aslMATFormat.cxx
	)

	add_library(aslmatio ${aslmatio_SOURCES})
	set_target_properties(aslmatio PROPERTIES VERSION ${ASL_VERSION} SOVERSION ${ASL_VERSION_MAJOR})
	target_link_libraries(aslmatio ${MATIO_LIBRARIES})
	INSTALL_SUBLIB(aslmatio aslmatio_PUBLIC_HEADERS)
endif (WITH_MATIO)

