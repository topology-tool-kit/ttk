cmake_minimum_required(VERSION 3.5)

project(contourForestsGui)

add_executable(contourForestsGui main.cpp)
target_link_libraries(contourForestsGui ttkContourForests ttkUserInterfaceBase)
install(TARGETS contourForestsGui RUNTIME DESTINATION bin/ttk)
