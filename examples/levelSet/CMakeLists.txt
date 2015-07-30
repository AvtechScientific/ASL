include_directories(${CMAKE_SOURCE_DIR}/src)


add_executable(levelSetBasic levelSetBasic.cc)
target_link_libraries(levelSetBasic aslnum aslvtk asl)
set_property(TARGET levelSetBasic PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/levelSetBasic)
INSTALL_EXAMPLE(levelSetBasic levelSetBasic.cc)

add_executable(levelSetNormalGrowth levelSetNormalGrowth.cc)
target_link_libraries(levelSetNormalGrowth aslnum aslvtk asl)
set_property(TARGET levelSetNormalGrowth PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/levelSetNormalGrowth)
INSTALL_EXAMPLE(levelSetNormalGrowth levelSetNormalGrowth.cc)

add_executable(levelSetFacetedGrowth levelSetFacetedGrowth.cc)
target_link_libraries(levelSetFacetedGrowth aslnum aslvtk asl)
set_property(TARGET levelSetFacetedGrowth PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/levelSetFacetedGrowth)
INSTALL_EXAMPLE(levelSetFacetedGrowth levelSetFacetedGrowth.cc)