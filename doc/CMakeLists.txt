if (WITH_INTERNAL_DOC)
	find_package(Doxygen REQUIRED)

	if (DOXYGEN_DOT_FOUND)
		set(DOXYGEN_HAVE_DOT "YES")
	else (DOXYGEN_DOT_FOUND)
		set(DOXYGEN_HAVE_DOT "NO")
		message(WARNING "Graphviz dot was not found! Class inheritance diagrams in the generated documentation will be limited.")
	endif (DOXYGEN_DOT_FOUND)

	# Configure Doxyfile
	configure_file(
		"${PROJECT_SOURCE_DIR}/doc/Developer-Guide/Doxyfile.in"
		"${PROJECT_BINARY_DIR}/doc/Developer-Guide/Doxyfile"
	)
	add_custom_target(docs ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/Developer-Guide/Doxyfile)
	install(CODE "execute_process(COMMAND ${CMAKE_BUILD_TOOL} docs)")
	install(DIRECTORY ${CMAKE_BINARY_DIR}/doc/html DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()