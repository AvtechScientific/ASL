include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(jumpingBox jumpingBox.cc)
target_link_libraries(jumpingBox aslnum aslvtk asl)
INSTALL_EXAMPLE(jumpingBox jumpingBox.cc)