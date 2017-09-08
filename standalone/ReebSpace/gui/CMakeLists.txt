cmake_minimum_required(VERSION 3.5)

project(reebSpaceGui)
add_executable(reebSpaceGui main.cpp)
target_link_libraries(reebSpaceGui ttkReebSpace ttkUserInterfaceBase)
install(TARGETS reebSpaceGui RUNTIME DESTINATION bin/ttk)

