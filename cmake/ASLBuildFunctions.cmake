#
# Advanced Simulation Library <http://asl.org.il>
# 
# Copyright 2015 Avtech Scientific <http://avtechscientific.com>
#
#
# This file is part of Advanced Simulation Library (ASL).
#
# ASL is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, version 3 of the License.
#
# ASL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with ASL. If not, see <http://www.gnu.org/licenses/>.
#


# Installs sublibraries - binaries and public headers
# (preserving source tree structure)
function(INSTALL_SUBLIB _SUBLIB _SUBLIB_PUBLIC_HEADERS)

	# Add current sublibrary to the list of all sublibs for inclusion in ASL.pc
	# using global property ASL_SUBLIBS_GLOBAL_PROPERTY
	set_property(GLOBAL APPEND_STRING PROPERTY ASL_SUBLIBS_GLOBAL_PROPERTY " -l${_SUBLIB}")

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
		INCLUDES DESTINATION ${ASL_INSTALL_INCLUDEDIR}
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
			${ASL_INSTALL_INCLUDEDIR}/${directories}
			COMPONENT Devel
		)
	endforeach()
endfunction(INSTALL_SUBLIB)


# Installs examples: binary and the corresponding source code (preserving source tree structure)
# puts `asl-` preffix before binaries in order to avoid name collisions
# creates a separate directory in the build tree to make experimenting easier
function(INSTALL_EXAMPLE _TARGET _SOURCE)
	# Add `asl-` preffix in order to avoid name collisions
	set_property(TARGET ${_TARGET} PROPERTY OUTPUT_NAME "asl-${_TARGET}")
	# Create a separate directory in the build tree for each example to make experimenting easier
	set_property(TARGET ${_TARGET} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_TARGET})
	install(TARGETS ${_TARGET} RUNTIME DESTINATION  ${CMAKE_INSTALL_BINDIR})

	# Determine relative path from ${CMAKE_SOURCE_DIR} to ${_SOURCE}
	file(RELATIVE_PATH relative_path ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${_SOURCE})
	# Extract directories of the relative path
	get_filename_component(directories ${relative_path} DIRECTORY)
	install(FILES ${_SOURCE} DESTINATION ${CMAKE_INSTALL_DOCDIR}/${directories})
endfunction(INSTALL_EXAMPLE)