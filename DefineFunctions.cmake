# Installs sublibraries - binaries and
# public headers (preserving source tree structure)
function(INSTALL_SUBLIB _SUBLIB _SUBLIB_PUBLIC_HEADERS)
	# ToDo: add COMPONENT headers
	install(TARGETS
		${_SUBLIB}
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib${LIB_SUFFIX}
		ARCHIVE DESTINATION lib${LIB_SUFFIX}
	)

	# Install public headers preserving the source tree structure
	foreach(header ${${_SUBLIB_PUBLIC_HEADERS}})
		# Determine relative path from ${CMAKE_SOURCE_DIR}/src to ${header}
		file(RELATIVE_PATH relative_path ${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/${header})
		# Extract directories of the relative path
		get_filename_component(directories ${relative_path} DIRECTORY)
		# ToDo: add COMPONENT libraries
		install(FILES
			${header}
			DESTINATION
			include/asl/${directories}
		)
	endforeach()
endfunction(INSTALL_SUBLIB)
