# ToDo: recheck here and everywhere whether all listed
# _PUBLIC_HEADERS are indeed needed and if not move them
# to the _SOURCES (or _PRIVATE_HEADERS?)

add_subdirectory(acl)
add_subdirectory(data)
add_subdirectory(math)
add_subdirectory(num)
add_subdirectory(numExtend)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/acl)


# aslcommon

set(aslcommon_PUBLIC_HEADERS
	utilities/aslParametersManager.h
	utilities/aslTimeStamp.h
	aslUtilities.h
)

set(aslcommon_SOURCES
	${aslcommon_PUBLIC_HEADERS}
	utilities/aslParametersManager.cxx
	utilities/aslTimeStamp.cxx
	aslUtilities.cxx
)

add_library(aslcommon ${aslcommon_SOURCES})
target_link_libraries(aslcommon PRIVATE ${Boost_LIBRARIES})
INSTALL_SUBLIB(aslcommon aslcommon_PUBLIC_HEADERS)


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
target_link_libraries(asl PUBLIC aslcommon aslacl asldata aslmath ${Boost_LIBRARIES})
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
target_link_libraries(aslvtk PRIVATE ${VTK_LIBRARIES} PRIVATE ${OpenCL_LIBRARIES})
INSTALL_SUBLIB(aslvtk aslvtk_PUBLIC_HEADERS)


# aslmatio

if (WITH_MATIO)
	find_package(MATIO 1.5.2 REQUIRED)
	include_directories(SYSTEM ${MATIO_INCLUDE_DIRS})

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
	target_link_libraries(aslmatio PRIVATE ${MATIO_LIBRARIES})
	INSTALL_SUBLIB(aslmatio aslmatio_PUBLIC_HEADERS)
endif (WITH_MATIO)

