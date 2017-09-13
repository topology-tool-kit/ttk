cmake_minimum_required(VERSION 3.2)

project(meshSubdivisionGui)
add_executable(meshSubdivisionGui main.cpp)
target_link_libraries(meshSubdivisionGui ttkMeshSubdivision ttkUserInterfaceBase)

install(TARGETS meshSubdivisionGui RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

