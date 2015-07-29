include_directories(${CMAKE_SOURCE_DIR}/src)

# Install asl-hardware
add_executable(asl-hardware asl-hardware.cc)
target_link_libraries(asl-hardware aslacl)
install(TARGETS asl-hardware RUNTIME DESTINATION  ${CMAKE_INSTALL_BINDIR})