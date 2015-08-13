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


# acloperators

set(acloperators_PUBLIC_HEADERS
	Operators/aclElementAssignmentSafe.h
	Operators/aclElementIfElse.h
	Operators/aclElementSum.h
	Operators/aclElementSubtraction.h
	Operators/aclElementSqrt.h
	Operators/aclElementSin.h
	Operators/aclElementProduct.h
	Operators/aclElementDivision.h
	Operators/aclElementGenericBinary.h
	Operators/aclGenericAtomicFunction.h
	Operators/aclElementGenericUnary.h
	Operators/aclElementCos.h
	Operators/aclElementSelect.h
	Operators/aclElementMad.h
	Operators/aclElementExcerpt.h
	Operators/aclElementParser.h
	Operators/aclElementSyncCopy.h
	Operators/aclOperatorUnary.h
	Operators/aclOperatorBinary.h
	Operators/aclOperatorTernary.h
	Operators/aclOperatorGeneric.h
	Operators/aclElementFor.h
	Operators/aclElementConvert.h
)

set(acloperators_SOURCES
	Operators/aclElementAssignmentSafe.cxx
	Operators/aclElementIfElse.cxx
	Operators/aclElementSum.cxx
	Operators/aclElementSubtraction.cxx
	Operators/aclElementSqrt.cxx
	Operators/aclElementSin.cxx
	Operators/aclElementProduct.cxx
	Operators/aclElementDivision.cxx
	Operators/aclElementGenericBinary.cxx
	Operators/aclGenericAtomicFunction.cxx
	Operators/aclElementGenericUnary.cxx
	Operators/aclElementCos.cxx
	Operators/aclElementSelect.cxx
	Operators/aclElementMad.cxx
	Operators/aclElementExcerpt.cxx
	Operators/aclElementParser.cxx
	Operators/aclElementSyncCopy.cxx
	Operators/aclOperatorUnary.cxx
	Operators/aclOperatorBinary.cxx
	Operators/aclOperatorTernary.cxx
	Operators/aclOperatorGeneric.cxx
	Operators/aclElementFor.cxx
	Operators/aclElementConvert.cxx
)
