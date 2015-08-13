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

add_executable(testSMDiff testSMDiff.cc)
target_link_libraries(testSMDiff aslnum aslvtk asl)
INSTALL_EXAMPLE(testSMDiff testSMDiff.cc)

add_executable(testSMDiff3C testSMDiff3C.cc)
target_link_libraries(testSMDiff3C aslnum aslvtk asl)
INSTALL_EXAMPLE(testSMDiff3C testSMDiff3C.cc)

add_executable(testSMPhi testSMPhi.cc)
target_link_libraries(testSMPhi aslnum aslvtk asl)
INSTALL_EXAMPLE(testSMPhi testSMPhi.cc)

add_executable(testSMPhiBV testSMPhiBV.cc)
target_link_libraries(testSMPhiBV aslnum aslvtk asl)
INSTALL_EXAMPLE(testSMPhiBV testSMPhiBV.cc)
