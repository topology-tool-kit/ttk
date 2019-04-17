cmake_minimum_required(VERSION 3.2)

project(ttkPointMergerGui)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(ttkPointMergerGui main.cpp)
target_link_libraries(ttkPointMergerGui ${VTK_LIBRARIES}
                      ttkPointMerger ttkProgramBase ttkUserInterfaceBase)

install(TARGETS ttkPointMergerGui RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})
