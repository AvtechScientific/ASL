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