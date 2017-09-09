cmake_minimum_required(VERSION 3.2)

project(geometrySmootherGui)

add_executable(geometrySmootherGui main.cpp)
target_link_libraries(geometrySmootherGui ttkGeometrySmoother ttkUserInterfaceBase)
install(TARGETS geometrySmootherGui RUNTIME DESTINATION bin/ttk)

