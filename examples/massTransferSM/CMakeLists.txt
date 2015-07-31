include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(testSMDiff testSMDiff.cc)
target_link_libraries(testSMDiff aslnum aslvtk asl)
INSTALL_EXAMPLE(testSMDiff testSMDiff.cc)

add_executable(testSMDiff3C testSMDiff3C.cc)
target_link_libraries(testSMDiff3C aslnum aslvtk asl)
INSTALL_EXAMPLE(testSMDiff3C testSMDiff3C.cc)

add_executable(testSMPhi testSMPhi.cc)
target_link_libraries(testSMPhi aslnum aslvtk asl)
INSTALL_EXAMPLE(testSMPhi testSMPhi.cc)

add_executable(testSMPhiBV testSMPhiBV.cc)
target_link_libraries(testSMPhiBV aslnum aslvtk asl)
INSTALL_EXAMPLE(testSMPhiBV testSMPhiBV.cc)