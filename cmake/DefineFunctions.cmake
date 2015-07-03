# Installs sublibraries - binaries and public headers
# (preserving source tree structure)
function(INSTALL_SUBLIB _SUBLIB _SUBLIB_PUBLIC_HEADERS)

	set_target_properties(
		${_SUBLIB} PROPERTIES VERSION ${ASL_VERSION}
		SOVERSION ${ASL_VERSION_MAJOR}
		INTERFACE_${_SUBLIB}_MAJOR_VERSION ${ASL_VERSION_MAJOR}
		COMPATIBLE_INTERFACE_STRING ${ASL_VERSION_MAJOR}
	)

	install(TARGETS
		${_SUBLIB} EXPORT ASLTargets
		RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
		INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/asl
	)

	# Install public headers preserving the source tree structure
	foreach(header ${${_SUBLIB_PUBLIC_HEADERS}})
		# Determine relative path from ${CMAKE_SOURCE_DIR}/src to ${header}
		file(RELATIVE_PATH relative_path ${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/${header})
		# Extract directories of the relative path
		get_filename_component(directories ${relative_path} DIRECTORY)
		install(FILES
			${header}
			DESTINATION
			${CMAKE_INSTALL_INCLUDEDIR}/asl/${directories}
			COMPONENT Devel
		)
	endforeach()
endfunction(INSTALL_SUBLIB)