include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(testVectorOfElements testVectorOfElements.cc)
target_link_libraries(testVectorOfElements aslacl)

add_executable(testMatrixOfElements testMatrixOfElements.cc)
target_link_libraries(testMatrixOfElements aslacl aslmath asldata asl)

add_executable(testPerformance testPerformance.cc)
target_link_libraries(testPerformance aslacl)

add_executable(testKernel testKernel.cc)
target_link_libraries(testKernel aslacl)

add_executable(testHardware testHardware.cc)
target_link_libraries(testHardware aslacl)

add_executable(testOperators testOperators.cc)
target_link_libraries(testOperators aslacl)

add_executable(testKernelMerger testKernelMerger.cc)
target_link_libraries(testKernelMerger aslacl)

add_executable(testPrivateVar testPrivateVar.cc)
target_link_libraries(testPrivateVar aslacl)
