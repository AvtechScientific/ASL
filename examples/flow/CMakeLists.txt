include_directories(${CMAKE_SOURCE_DIR}/src)


add_executable(flow flow.cc)
target_link_libraries(flow aslnum aslvtk asl)
set_property(TARGET flow PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flow)

add_executable(flow2 flow2.cc)
target_link_libraries(flow2 aslnum aslvtk asl)
set_property(TARGET flow2 PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flow2)

add_executable(flow3 flow3.cc)
target_link_libraries(flow3 aslnum aslvtk asl)
set_property(TARGET flow3 PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flow3)

add_executable(flowRotatingCylinders flowRotatingCylinders.cc)
target_link_libraries(flowRotatingCylinders aslnum aslvtk asl)
set_property(TARGET flowRotatingCylinders PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flowRotatingCylinders)

add_executable(flowKDPGrowth flowKDPGrowth.cc)
target_link_libraries(flowKDPGrowth aslnum aslvtk asl)
set_property(TARGET flowKDPGrowth PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flowKDPGrowth)

add_executable(locomotive_laminar locomotive_laminar.cc)
target_link_libraries(locomotive_laminar aslnum aslvtk asl)
set_property(TARGET locomotive_laminar PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/locomotive_laminar)

add_executable(compressor compressor.cc)
target_link_libraries(compressor aslnum aslvtk asl)
set_property(TARGET compressor PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/compressor)

add_executable(locomotive_in_tunnel locomotive_in_tunnel.cc)
target_link_libraries(locomotive_in_tunnel aslnum aslvtk asl)
set_property(TARGET locomotive_in_tunnel PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/locomotive_in_tunnel)

add_executable(locomotive_stability locomotive_stability.cc)
target_link_libraries(locomotive_stability aslnum aslvtk asl)
set_property(TARGET locomotive_stability PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/locomotive_stability)

add_executable(bus_wind bus_wind.cc)
target_link_libraries(bus_wind aslnum aslvtk asl)
set_property(TARGET bus_wind PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bus_wind)

add_executable(multiphase_flow multiphase_flow.cc)
target_link_libraries(multiphase_flow aslnum aslvtk asl)
set_property(TARGET multiphase_flow PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/multiphase_flow)

add_executable(multicomponent_flow multicomponent_flow.cc)
target_link_libraries(multicomponent_flow aslnum aslvtk asl)
set_property(TARGET multicomponent_flow PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/multicomponent_flow)

add_executable(pitot_tube_ice pitot_tube_ice.cc)
target_link_libraries(pitot_tube_ice aslnum aslvtk asl)
set_property(TARGET pitot_tube_ice PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pitot_tube_ice)
