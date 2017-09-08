cmake_minimum_required(VERSION 3.5)

project(contourForestsCmd)

add_executable(contourForestsCmd main.cpp)
target_link_libraries(contourForestsCmd ttkContourForests ttkProgramBase)
install(TARGETS contourForestsCmd RUNTIME DESTINATION bin/ttk)
