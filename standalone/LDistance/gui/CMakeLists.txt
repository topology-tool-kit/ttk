cmake_minimum_required(VERSION 3.5)

project(lDistanceGui)
add_executable(lDistanceGui main.cpp)
target_link_libraries(lDistanceGui ttkLDistance ttkUserInterfaceBase)

