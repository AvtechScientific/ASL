include_directories(${CMAKE_SOURCE_DIR}/src)


add_executable(poroelastic poroelastic.cc)
target_link_libraries(poroelastic asl aslnum aslvtk)
set_property(TARGET poroelastic PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/poroelastic)

add_executable(cubeIncompressibleGravity cubeIncompressibleGravity.cc)
target_link_libraries(cubeIncompressibleGravity asl aslnum aslvtk)
set_property(TARGET cubeIncompressibleGravity PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cubeIncompressibleGravity)

add_executable(cubePoroelasticGravity cubePoroelasticGravity.cc)
target_link_libraries(cubePoroelasticGravity asl aslnum aslvtk)
set_property(TARGET cubePoroelasticGravity PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cubePoroelasticGravity)

add_executable(cubeGravity cubeGravity.cc)
target_link_libraries(cubeGravity asl aslnum aslvtk)
set_property(TARGET cubeGravity PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cubeGravity)

add_executable(acousticWaves acousticWaves.cc)
target_link_libraries(acousticWaves asl aslnum aslvtk)
set_property(TARGET acousticWaves PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/acousticWaves)

