cmake_minimum_required(VERSION 3.5)

project(geometrySmootherCmd)

add_executable(geometrySmootherCmd main.cpp)
target_link_libraries(geometrySmootherCmd scalarFieldSmoother ttkProgramBase)
install(TARGETS geometrySmootherCmd RUNTIME DESTINATION bin/ttk)

