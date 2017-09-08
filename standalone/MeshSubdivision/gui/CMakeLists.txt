cmake_minimum_required(VERSION 3.5)

project(meshSubdivisionGui)
add_executable(meshSubdivisionGui main.cpp)
target_link_libraries(meshSubdivisionGui ttkMeshSubdivision ttkUserInterfaceBase)
install(TARGETS meshSubdivisionGui RUNTIME DESTINATION bin/ttk)

