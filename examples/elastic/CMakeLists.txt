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


add_executable(poroelastic poroelastic.cc)
target_link_libraries(poroelastic aslnum aslvtk asl)
INSTALL_EXAMPLE(poroelastic poroelastic.cc)

add_executable(cubeIncompressibleGravity cubeIncompressibleGravity.cc)
target_link_libraries(cubeIncompressibleGravity aslnum aslvtk asl)
INSTALL_EXAMPLE(cubeIncompressibleGravity cubeIncompressibleGravity.cc)

add_executable(cubePoroelasticGravity cubePoroelasticGravity.cc)
target_link_libraries(cubePoroelasticGravity aslnum aslvtk asl)
INSTALL_EXAMPLE(cubePoroelasticGravity cubePoroelasticGravity.cc)

add_executable(cubeGravity cubeGravity.cc)
target_link_libraries(cubeGravity aslnum aslvtk asl)
INSTALL_EXAMPLE(cubeGravity cubeGravity.cc)

add_executable(acousticWaves acousticWaves.cc)
target_link_libraries(acousticWaves aslnum aslvtk asl)
INSTALL_EXAMPLE(acousticWaves acousticWaves.cc)
