cmake_minimum_required(VERSION 3.5)

project(lDistanceCmd)

add_executable(lDistanceCmd main.cpp)
target_link_libraries(lDistanceCmd ttkLDistance ttkProgramBase)
install(TARGETS lDistanceCmd RUNTIME DESTINATION bin/ttk)

