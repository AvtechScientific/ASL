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


# aclmath

set(aclmath_PUBLIC_HEADERS
	aclMath/aclVectorOfElements.h
	aclMath/aclVectorOfElementsDef.h
	aclMath/aclVectorOfElementsOperations.h
	aclMath/aclMatrixOfElements.h
	aclMath/aclComplexNumOfElements.h
	aclMath/aclQuaternionOfElements.h
	aclMath/aclReductionAlgGenerator.h
	aclMath/aclBarycentric.h
	aclMath/aclMathAlg.h
)

set(aclmath_SOURCES
	aclMath/aclVectorOfElementsDef.cxx
	aclMath/aclVectorOfElementsOperations.cxx
	aclMath/aclMatrixOfElements.cxx
	aclMath/aclComplexNumOfElements.cxx
	aclMath/aclQuaternionOfElements.cxx
	aclMath/aclReductionAlgGenerator.cxx
	aclMath/aclBarycentric.cxx
	aclMath/aclMathAlg.cxx
)
