cmake_minimum_required(VERSION 3.2)

project(reebSpaceCmd)
add_executable(reebSpaceCmd main.cpp)
target_link_libraries(reebSpaceCmd ttkReebSpace ttkProgramBase)
install(TARGETS reebSpaceCmd RUNTIME DESTINATION bin/ttk)

