
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(testASLData testASLData.cc)
target_link_libraries(testASLData aslmath asl)

add_executable(testDistanceFunction testDistanceFunction.cc)
target_link_libraries(testDistanceFunction aslmath aslvtk aslnum asl)

add_executable(testReductionFunction testReductionFunction.cc)
target_link_libraries(testReductionFunction asl)
