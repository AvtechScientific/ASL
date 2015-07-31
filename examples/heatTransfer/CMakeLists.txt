include_directories(${CMAKE_SOURCE_DIR}/src)


add_executable(surfaceFlux surfaceFlux.cc)
target_link_libraries(surfaceFlux aslnum aslvtk asl)
INSTALL_EXAMPLE(surfaceFlux surfaceFlux.cc)