cmake_minimum_required(VERSION 3.5)

# name of the project
project(ttkExample-c++)

set(CMAKE_CXX_STANDARD 11)

find_package(TTKBase REQUIRED)

add_executable(ttkExample-c++ main.cpp)

target_link_libraries(ttkExample-c++
  PUBLIC
    ttk::base::baseAll
    )
