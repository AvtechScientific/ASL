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


# acldatatypes

set(acldatatypes_PUBLIC_HEADERS
	DataTypes/aclConstant.h
	DataTypes/aclIndex.h 
	DataTypes/aclGroupID.h
	DataTypes/aclIndexExt.h 
	DataTypes/aclPrivateVariable.h 
	DataTypes/aclPrivateArray.h
	DataTypes/aclLocalArray.h
	DataTypes/aclVariable.h 
	DataTypes/aclVariableReference.h 
	DataTypes/aclVariableSP.h
	DataTypes/aclMemBlock.h
	DataTypes/aclArray.h
	DataTypes/aclSubvector.h
)

set(acldatatypes_SOURCES
	DataTypes/aclConstant.cxx
	DataTypes/aclIndex.cxx
	DataTypes/aclGroupID.cxx
	DataTypes/aclIndexExt.cxx
	DataTypes/aclPrivateVariable.cxx 
	DataTypes/aclPrivateArray.cxx
	DataTypes/aclLocalArray.cxx
	DataTypes/aclVariable.cxx 
	DataTypes/aclVariableReference.cxx 
	DataTypes/aclVariableSP.cxx 
	DataTypes/aclMemBlock.cxx
	DataTypes/aclArray.cxx
	DataTypes/aclSubvector.cxx
)
