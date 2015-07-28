
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(testVTK_IO testVTK_IO.cc)
target_link_libraries(testVTK_IO aslnum aslvtk asl)
add_test(NAME testVTK_IO COMMAND testVTK_IO)
set_tests_properties(testVTK_IO PROPERTIES LABELS IO)