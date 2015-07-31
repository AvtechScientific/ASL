add_subdirectory(flow)
add_subdirectory(elastic)
add_subdirectory(levelSet)
add_subdirectory(jumpingObjects)
add_subdirectory(heatTransfer)
add_subdirectory(massTransferSM)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/examples/input_data DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples)