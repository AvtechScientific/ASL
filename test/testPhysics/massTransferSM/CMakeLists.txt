include_directories(${CMAKE_SOURCE_DIR}/src)


add_executable(testSMDiff testSMDiff.cc)
target_link_libraries(testSMDiff asl aslnum aslvtk)
set_property(TARGET testSMDiff PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testSMDiff)

add_executable(testSMDiff3C testSMDiff3C.cc)
target_link_libraries(testSMDiff3C asl aslnum aslvtk)
set_property(TARGET testSMDiff3C PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testSMDiff3C)

add_executable(testSMPhi testSMPhi.cc)
target_link_libraries(testSMPhi asl aslnum aslvtk)
set_property(TARGET testSMPhi PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testSMPhi)

add_executable(testSMPhiBV testSMPhiBV.cc)
target_link_libraries(testSMPhiBV asl aslnum aslvtk)
set_property(TARGET testSMPhiBV PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testSMPhiBV)

