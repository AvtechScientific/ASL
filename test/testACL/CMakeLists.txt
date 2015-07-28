include_directories(${CMAKE_SOURCE_DIR}/src)

# Install asl-hardware
add_executable(asl-hardware asl-hardware.cc)
target_link_libraries(asl-hardware aslacl)
install(TARGETS asl-hardware RUNTIME DESTINATION  ${CMAKE_INSTALL_BINDIR})

add_executable(testVectorOfElements testVectorOfElements.cc)
target_link_libraries(testVectorOfElements aslacl)
add_test(NAME testVectorOfElements COMMAND testVectorOfElements)

add_executable(testMatrixOfElements testMatrixOfElements.cc)
target_link_libraries(testMatrixOfElements aslacl aslmath asldata asl)
add_test(NAME testMatrixOfElements COMMAND testMatrixOfElements)

add_executable(testPerformance testPerformance.cc)
target_link_libraries(testPerformance aslacl)
add_test(NAME testPerformance COMMAND testPerformance)
set_tests_properties(testPerformance PROPERTIES LABELS Performance)

add_executable(testKernel testKernel.cc)
target_link_libraries(testKernel aslacl)
add_test(NAME testKernel COMMAND testKernel)

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