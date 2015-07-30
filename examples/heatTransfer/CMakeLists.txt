include_directories(${CMAKE_SOURCE_DIR}/src)


add_executable(surfaceFlux surfaceFlux.cc)
target_link_libraries(surfaceFlux aslnum aslvtk asl)
set_property(TARGET surfaceFlux PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/surfaceFlux)
INSTALL_EXAMPLE(surfaceFlux surfaceFlux.cc)