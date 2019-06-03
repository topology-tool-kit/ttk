cmake_minimum_required(VERSION 3.8)

# name of the project
project(ttkExample-vtk-c++ LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)

find_package(TTKVTK REQUIRED)

add_executable(ttkExample-vtk-c++ main.cpp)

target_link_libraries(ttkExample-vtk-c++
  PUBLIC
    ttk::vtk::ttkAll
    )
