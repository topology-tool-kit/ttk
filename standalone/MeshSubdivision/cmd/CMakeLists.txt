cmake_minimum_required(VERSION 3.2)

project(meshSubdivisionCmd)
add_executable(meshSubdivisionCmd main.cpp)
target_link_libraries(meshSubdivisionCmd ttkMeshSubdivision ttkProgramBase)

install(TARGETS meshSubdivisionCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

