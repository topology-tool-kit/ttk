cmake_minimum_required(VERSION 3.5)
project(jacobiSetCmd)

add_executable(jacobiSetCmd main.cpp)
target_link_libraries(jacobiSetCmd ttkJacobiSet ttkProgramBase)

