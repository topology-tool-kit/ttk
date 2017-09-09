cmake_minimum_required(VERSION 3.2)

project(contourForestsCmd)

add_executable(contourForestsCmd main.cpp)
target_link_libraries(contourForestsCmd ttkContourForests ttkProgramBase)
install(TARGETS contourForestsCmd RUNTIME DESTINATION bin/ttk)
