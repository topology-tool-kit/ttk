cmake_minimum_required(VERSION 3.2)

project(jacobiSetGui)

add_executable(jacobiSetGui main.cpp)
target_link_libraries(jacobiSetGui ttkJacobiSet ttkUserInterfaceBase)
install(TARGETS jacobiSetGui RUNTIME DESTINATION bin/ttk)

