include_directories(${CMAKE_SOURCE_DIR}/src)

# Install asl-hardware
add_executable(asl-hardware asl-hardware.cc)
target_link_libraries(asl-hardware aslacl)
install(TARGETS asl-hardware RUNTIME DESTINATION  ${CMAKE_INSTALL_BINDIR})

add_executable(testVectorOfElements testVectorOfElements.cc)
target_link_libraries(testVectorOfElements aslacl)
add_test(testVectorOfElements, testVectorOfElements)

add_executable(testMatrixOfElements testMatrixOfElements.cc)
target_link_libraries(testMatrixOfElements aslacl aslmath asldata asl)
add_test(testMatrixOfElements, testMatrixOfElements)

add_executable(testPerformance testPerformance.cc)
target_link_libraries(testPerformance aslacl)
add_test(testPerformance, testPerformance)

add_executable(testKernel testKernel.cc)
target_link_libraries(testKernel aslacl)
add_test(testKernel, testKernel)

add_executable(testOperators testOperators.cc)
target_link_libraries(testOperators aslacl)
add_test(testOperators, testOperators)

add_executable(testKernelMerger testKernelMerger.cc)
target_link_libraries(testKernelMerger aslacl)
add_test(testKernelMerger, testKernelMerger)

add_executable(testPrivateVar testPrivateVar.cc)
target_link_libraries(testPrivateVar aslacl)
add_test(testPrivateVar, testPrivateVar)