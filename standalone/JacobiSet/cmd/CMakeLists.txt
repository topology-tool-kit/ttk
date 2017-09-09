cmake_minimum_required(VERSION 3.2)
project(jacobiSetCmd)

add_executable(jacobiSetCmd main.cpp)
target_link_libraries(jacobiSetCmd ttkJacobiSet ttkProgramBase)
install(TARGETS jacobiSetCmd RUNTIME DESTINATION bin/ttk)

