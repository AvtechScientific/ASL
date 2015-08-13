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


include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/acl)


# aslnum

set(aslnum_PUBLIC_HEADERS
	aslNumMethod.h
	aslSingleKernelNM.h	
	aslFDAdvectionDiffusion.h
	aslFDAdvectionDiffusionBC.h
	aslFDStefanMaxwell.h
	aslFDElChemKinetics.h
	aslFDMultiPhase.h
	aslFDAdvectionDiffusion2.h
	aslFDAdvectionDiffusionInhomogeneous.h
	aslLBGK.h
	aslLBGKBC.h
	aslFDElasticity.h
	aslFDPoroElasticity.h
	aslFDPoroElasticityBC.h
	aslBCond.h
	aslBasicBC.h
	aslBasicBC2.h
	aslCrystalGrowthBC.h
	aslNumMethodsMerger.h
	aslFDElasticityBC.h
	aslLevelSet.h
	aslLevelSetLinear.h
	aslDFOptimizer.h
	aslLSNormalGrowth.h
	aslInterfaceTrackingAlg1.h
	aslLSFacetedGrowth.h
	aslTimeContinuations.h
	aslDataResampling.h
)

set(aslnum_SOURCES
	aslNumMethod.cxx
	aslSingleKernelNM.cxx
	aslFDAdvectionDiffusion.cxx
	aslFDAdvectionDiffusionBC.cxx
	aslFDStefanMaxwell.cxx
	aslFDElChemKinetics.cxx
	aslFDMultiPhase.cxx
	aslFDAdvectionDiffusion2.cxx
	aslFDAdvectionDiffusionInhomogeneous.cxx
	aslLBGK.cxx
	aslLBGKBC.cxx
	aslFDElasticity.cxx
	aslFDPoroElasticity.cxx
	aslFDPoroElasticityBC.cxx
	aslBCond.cxx
	aslBasicBC.cxx
	aslBasicBC2.cxx
	aslCrystalGrowthBC.cxx
	aslFDElasticityBC.cxx
	aslNumMethodsMerger.cxx
	aslLevelSet.cxx
	aslLevelSetLinear.cxx
	aslDFOptimizer.cxx
	aslInterfaceTrackingAlg1.cxx
	aslLSNormalGrowth.cxx
	aslLSFacetedGrowth.cxx
	aslTimeContinuations.cxx
	aslDataResampling.cxx
)

add_library(aslnum ${aslnum_PUBLIC_HEADERS} ${aslnum_SOURCES})
target_link_libraries(aslnum PRIVATE asl aslacl aslcommon ${OpenCL_LIBRARIES})
INSTALL_SUBLIB(aslnum aslnum_PUBLIC_HEADERS)
