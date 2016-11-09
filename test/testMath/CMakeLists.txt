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


include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(testASLData testASLData.cc)
target_link_libraries(testASLData asl)
add_test(NAME testASLData COMMAND testASLData)

add_executable(testDistanceFunction testDistanceFunction.cc)
target_link_libraries(testDistanceFunction aslvtk aslnum asl)
add_test(NAME testDistanceFunction COMMAND testDistanceFunction)
set_tests_properties(testDistanceFunction PROPERTIES LABELS DoublePrecision)

add_executable(testReductionFunction testReductionFunction.cc)
target_link_libraries(testReductionFunction asl)
add_test(NAME testReductionFunction COMMAND testReductionFunction)
set_tests_properties(testReductionFunction PROPERTIES LABELS DoublePrecision)