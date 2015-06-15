if (WITH_INTERNAL_DOC)
	find_package(Doxygen REQUIRED)
	# Configure Doxyfile
	configure_file(
		"${PROJECT_SOURCE_DIR}/doc/Developer-Guide/Doxyfile.in"
		"${PROJECT_BINARY_DIR}/doc/Developer-Guide/Doxyfile"
	)
	add_custom_target(docs ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/Developer-Guide/Doxyfile)
	install(CODE "execute_process(COMMAND ${CMAKE_BUILD_TOOL} docs)")
	install(DIRECTORY ${CMAKE_BINARY_DIR}/doc/html DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()
