include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(testSMDiff testSMDiff.cc)
target_link_libraries(testSMDiff aslnum aslvtk asl)
set_property(TARGET testSMDiff PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testSMDiff)
INSTALL_EXAMPLE(testSMDiff testSMDiff.cc)

add_executable(testSMDiff3C testSMDiff3C.cc)
target_link_libraries(testSMDiff3C aslnum aslvtk asl)
set_property(TARGET testSMDiff3C PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testSMDiff3C)
INSTALL_EXAMPLE(testSMDiff3C testSMDiff3C.cc)

add_executable(testSMPhi testSMPhi.cc)
target_link_libraries(testSMPhi aslnum aslvtk asl)
set_property(TARGET testSMPhi PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testSMPhi)
INSTALL_EXAMPLE(testSMPhi testSMPhi.cc)

add_executable(testSMPhiBV testSMPhiBV.cc)
target_link_libraries(testSMPhiBV aslnum aslvtk asl)
set_property(TARGET testSMPhiBV PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testSMPhiBV)
INSTALL_EXAMPLE(testSMPhiBV testSMPhiBV.cc)