
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(testVTK_IO testVTK_IO.cc)
target_link_libraries(testVTK_IO asl aslnum aslvtk)

