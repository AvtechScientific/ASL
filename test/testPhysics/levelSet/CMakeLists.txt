include_directories(${CMAKE_SOURCE_DIR}/src)


add_executable(levelSetBasic levelSetBasic.cc)
target_link_libraries(levelSetBasic asl aslnum aslvtk)
set_property(TARGET levelSetBasic PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/levelSetBasic)

add_executable(levelSetNormalGrowth levelSetNormalGrowth.cc)
target_link_libraries(levelSetNormalGrowth asl aslnum aslvtk)
set_property(TARGET levelSetNormalGrowth PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/levelSetNormalGrowth)

add_executable(levelSetFacetedGrowth levelSetFacetedGrowth.cc)
target_link_libraries(levelSetFacetedGrowth asl aslnum aslvtk)
set_property(TARGET levelSetFacetedGrowth PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/levelSetFacetedGrowth)

