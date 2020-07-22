cmake_minimum_required(VERSION 3.5)

project(bettiNumbers)

find_package(TTKVTK REQUIRED)

add_executable(betti main.cpp)

target_link_libraries(betti
  PUBLIC
    ttk::vtk::ttkAll
    )