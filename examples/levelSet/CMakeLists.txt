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


add_executable(levelSetBasic levelSetBasic.cc)
target_link_libraries(levelSetBasic aslnum aslvtk asl)
INSTALL_EXAMPLE(levelSetBasic levelSetBasic.cc)

add_executable(levelSetNormalGrowth levelSetNormalGrowth.cc)
target_link_libraries(levelSetNormalGrowth aslnum aslvtk asl)
INSTALL_EXAMPLE(levelSetNormalGrowth levelSetNormalGrowth.cc)

add_executable(levelSetFacetedGrowth levelSetFacetedGrowth.cc)
target_link_libraries(levelSetFacetedGrowth aslnum aslvtk asl)
INSTALL_EXAMPLE(levelSetFacetedGrowth levelSetFacetedGrowth.cc parameters_levelSetFacetedGrowth.ini)
