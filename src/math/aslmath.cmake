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


# aslmath

set(aslmath_PUBLIC_HEADERS
	math/aslVectors.h
	math/aslVectorsDynamicLength.h
	math/aslVectorsDynamicLengthOperations.h
	math/aslMatrices.h
	math/aslBarycentric.h
	math/aslInterpolation.h
	math/aslTemplates.h
	math/aslTemplatesExtras.h
	math/aslProbeTemplates.h
	math/aslTemplateVE.h
	math/aslTemplateVEExtras.h
	math/aslIndex2Position.h
	math/aslDistanceFunction.h
	math/aslDistanceFunctionAlg.h
	math/aslPositionFunction.h
)

set(aslmath_SOURCES
	math/aslMatrices.cxx
	math/aslBarycentric.cxx
	math/aslTemplates.cxx
	math/aslTemplatesExtras.cxx
	math/aslProbeTemplates.cxx
	math/aslTemplateVE.cxx
	math/aslTemplateVEExtras.cxx
	math/aslIndex2Position.cxx
	math/aslDistanceFunction.cxx
	math/aslDistanceFunctionAlg.cxx
	math/aslPositionFunction.cxx
)
