cmake_minimum_required(VERSION 3.2)

project(persistenceDiagramGui)
add_executable(persistenceDiagramGui main.cpp)
target_link_libraries(persistenceDiagramGui ttkPersistenceDiagram ttkUserInterfaceBase)
install(TARGETS persistenceDiagramGui RUNTIME DESTINATION bin/ttk)

