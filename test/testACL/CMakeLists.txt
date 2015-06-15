include_directories(${CMAKE_SOURCE_DIR}/src)

# Install asl-hardware
add_executable(asl-hardware testHardware.cc)
target_link_libraries(asl-hardware aslacl)
install(TARGETS asl-hardware RUNTIME DESTINATION  ${CMAKE_INSTALL_BINDIR})

add_executable(testVectorOfElements testVectorOfElements.cc)
target_link_libraries(testVectorOfElements aslacl)

add_executable(testMatrixOfElements testMatrixOfElements.cc)
target_link_libraries(testMatrixOfElements aslacl aslmath asldata asl)

add_executable(testPerformance testPerformance.cc)
target_link_libraries(testPerformance aslacl)

add_executable(testKernel testKernel.cc)
target_link_libraries(testKernel aslacl)

add_executable(testOperators testOperators.cc)
target_link_libraries(testOperators aslacl)

add_executable(testKernelMerger testKernelMerger.cc)
target_link_libraries(testKernelMerger aslacl)

add_executable(testPrivateVar testPrivateVar.cc)
target_link_libraries(testPrivateVar aslacl)
