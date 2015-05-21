include_directories(${CMAKE_SOURCE_DIR}/src)


add_executable(jumpingBox jumpingBox.cc)
target_link_libraries(jumpingBox asl aslnum aslvtk)
set_property(TARGET jumpingBox PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/jumpingBox)

