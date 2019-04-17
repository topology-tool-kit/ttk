cmake_minimum_required(VERSION 3.2)

project(ttkBlankCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(ttkBlankCmd main.cpp)
target_link_libraries(ttkBlankCmd ${VTK_LIBRARIES} ttkBlank ttkProgramBase)

install(TARGETS ttkBlankCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})
