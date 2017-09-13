cmake_minimum_required(VERSION 3.2)

project(scalarFieldSmootherCmd)
add_executable(scalarFieldSmootherCmd main.cpp)
target_link_libraries(scalarFieldSmootherCmd ttkScalarFieldSmoother ttkProgramBase)

install(TARGETS scalarFieldSmootherCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

