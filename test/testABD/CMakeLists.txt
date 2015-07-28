
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(testABDFormat testABDFormat.cc)
target_link_libraries(testABDFormat asldata aslmath asl)
add_test(testABDFormat, testABDFormat)