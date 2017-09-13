cmake_minimum_required(VERSION 3.2)

project(geometrySmootherCmd)

add_executable(geometrySmootherCmd main.cpp)
target_link_libraries(geometrySmootherCmd scalarFieldSmoother ttkProgramBase)

install(TARGETS geometrySmootherCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

