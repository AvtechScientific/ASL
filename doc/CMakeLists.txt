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
