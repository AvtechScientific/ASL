include_directories(${CMAKE_SOURCE_DIR}/src)


add_executable(levelSetBasic levelSetBasic.cc)
target_link_libraries(levelSetBasic aslnum aslvtk asl)
INSTALL_EXAMPLE(levelSetBasic levelSetBasic.cc)

add_executable(levelSetNormalGrowth levelSetNormalGrowth.cc)
target_link_libraries(levelSetNormalGrowth aslnum aslvtk asl)
INSTALL_EXAMPLE(levelSetNormalGrowth levelSetNormalGrowth.cc)

add_executable(levelSetFacetedGrowth levelSetFacetedGrowth.cc)
target_link_libraries(levelSetFacetedGrowth aslnum aslvtk asl)
INSTALL_EXAMPLE(levelSetFacetedGrowth levelSetFacetedGrowth.cc)