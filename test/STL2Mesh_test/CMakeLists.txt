set(cube_stl_ascii ${CMAKE_SOURCE_DIR}/test/STL2Mesh_test/cube_ascii.stl)
set(cube_stl_binary ${CMAKE_SOURCE_DIR}/test/STL2Mesh_test/cube_binary.stl)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/STL2Mesh_test_config.h.cmake ${CMAKE_BINARY_DIR}/test_inc/STL2Mesh_test_config.h )
include_directories(${CMAKE_BINARY_DIR}/test_inc)


add_smesh_test(STL2Mesh_test "SMESH")
