cmake_minimum_required(VERSION 3.2)

project(persistenceDiagramCmd)
add_executable(persistenceDiagramCmd main.cpp)
target_link_libraries(persistenceDiagramCmd ttkPersistenceDiagram ttkProgramBase)

install(TARGETS persistenceDiagramCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

