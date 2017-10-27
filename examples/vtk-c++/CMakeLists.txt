cmake_minimum_required(VERSION 3.5)

# name of the project
project(ttkExample-vtk-c++)

set(CMAKE_CXX_STANDARD 11)

find_package(VTK REQUIRED)
find_package(TTKBase REQUIRED)
find_package(TTKVTK REQUIRED)

include(${VTK_USE_FILE})

add_executable(ttkExample-vtk-c++ main.cpp)

target_link_libraries(ttkExample-vtk-c++ PUBLIC
  ${VTK_LIBRARIES}
  ttk::vtk::ttkAll)
