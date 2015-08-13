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


add_executable(flow flow.cc)
target_link_libraries(flow aslnum aslvtk asl)
INSTALL_EXAMPLE(flow flow.cc)

add_executable(flow2 flow2.cc)
target_link_libraries(flow2 aslnum aslvtk asl)
INSTALL_EXAMPLE(flow2 flow2.cc)

add_executable(flow3 flow3.cc)
target_link_libraries(flow3 aslnum aslvtk asl)
INSTALL_EXAMPLE(flow3 flow3.cc)

add_executable(flowRotatingCylinders flowRotatingCylinders.cc)
target_link_libraries(flowRotatingCylinders aslnum aslvtk asl)
INSTALL_EXAMPLE(flowRotatingCylinders flowRotatingCylinders.cc)

add_executable(flowKDPGrowth flowKDPGrowth.cc)
target_link_libraries(flowKDPGrowth aslnum aslvtk asl)
INSTALL_EXAMPLE(flowKDPGrowth flowKDPGrowth.cc)

add_executable(locomotive_laminar locomotive_laminar.cc)
target_link_libraries(locomotive_laminar aslnum aslvtk asl)
INSTALL_EXAMPLE(locomotive_laminar locomotive_laminar.cc)

add_executable(compressor compressor.cc)
target_link_libraries(compressor aslnum aslvtk asl)
INSTALL_EXAMPLE(compressor compressor.cc)

add_executable(locomotive locomotive.cc)
target_link_libraries(locomotive aslnum aslvtk asl)
INSTALL_EXAMPLE(locomotive locomotive.cc)

add_executable(locomotive_stability locomotive_stability.cc)
target_link_libraries(locomotive_stability aslnum aslvtk asl)
INSTALL_EXAMPLE(locomotive_stability locomotive_stability.cc)

add_executable(bus_wind bus_wind.cc)
target_link_libraries(bus_wind aslnum aslvtk asl)
INSTALL_EXAMPLE(bus_wind bus_wind.cc)

add_executable(multiphase_flow multiphase_flow.cc)
target_link_libraries(multiphase_flow aslnum aslvtk asl)
INSTALL_EXAMPLE(multiphase_flow multiphase_flow.cc)

add_executable(multicomponent_flow multicomponent_flow.cc)
target_link_libraries(multicomponent_flow aslnum aslvtk asl)
INSTALL_EXAMPLE(multicomponent_flow multicomponent_flow.cc)

add_executable(pitot_tube_ice pitot_tube_ice.cc)
target_link_libraries(pitot_tube_ice aslnum aslvtk asl)
INSTALL_EXAMPLE(pitot_tube_ice pitot_tube_ice.cc)
