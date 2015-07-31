include_directories(${CMAKE_SOURCE_DIR}/src)


add_executable(poroelastic poroelastic.cc)
target_link_libraries(poroelastic aslnum aslvtk asl)
INSTALL_EXAMPLE(poroelastic poroelastic.cc)

add_executable(cubeIncompressibleGravity cubeIncompressibleGravity.cc)
target_link_libraries(cubeIncompressibleGravity aslnum aslvtk asl)
INSTALL_EXAMPLE(cubeIncompressibleGravity cubeIncompressibleGravity.cc)

add_executable(cubePoroelasticGravity cubePoroelasticGravity.cc)
target_link_libraries(cubePoroelasticGravity aslnum aslvtk asl)
INSTALL_EXAMPLE(cubePoroelasticGravity cubePoroelasticGravity.cc)

add_executable(cubeGravity cubeGravity.cc)
target_link_libraries(cubeGravity aslnum aslvtk asl)
INSTALL_EXAMPLE(cubeGravity cubeGravity.cc)

add_executable(acousticWaves acousticWaves.cc)
target_link_libraries(acousticWaves aslnum aslvtk asl)
INSTALL_EXAMPLE(acousticWaves acousticWaves.cc)