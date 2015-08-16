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

add_executable(testVectorOfElements testVectorOfElements.cc)
target_link_libraries(testVectorOfElements aslacl)
add_test(NAME testVectorOfElements COMMAND testVectorOfElements)

add_executable(testMatrixOfElements testMatrixOfElements.cc)
target_link_libraries(testMatrixOfElements asl aslacl aslcommon)
add_test(NAME testMatrixOfElements COMMAND testMatrixOfElements)

# This test needs to be updated.
# add_executable(testPerformance testPerformance.cc)
# target_link_libraries(testPerformance aslacl)
# add_test(NAME testPerformance COMMAND testPerformance)
# set_tests_properties(testPerformance PROPERTIES LABELS Performance)

add_executable(testKernel testKernel.cc)
target_link_libraries(testKernel aslacl)
add_test(NAME testKernel COMMAND testKernel)
set_tests_properties(testKernel PROPERTIES LABELS DoublePrecision)

add_executable(testOperators testOperators.cc)
target_link_libraries(testOperators aslacl)
add_test(NAME testOperators COMMAND testOperators)

add_executable(testKernelMerger testKernelMerger.cc)
target_link_libraries(testKernelMerger aslacl)
add_test(NAME testKernelMerger COMMAND testKernelMerger)

add_executable(testPrivateVar testPrivateVar.cc)
target_link_libraries(testPrivateVar aslacl)
add_test(NAME testPrivateVar COMMAND testPrivateVar)
set_tests_properties(testPrivateVar PROPERTIES LABELS Performance)
