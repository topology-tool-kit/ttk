cmake_minimum_required(VERSION 3.2)

project(ttkManifoldCheckGui)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(ttkManifoldCheckGui main.cpp)
target_link_libraries(ttkManifoldCheckGui ${VTK_LIBRARIES}
                      ttkManifoldCheck ttkProgramBase ttkUserInterfaceBase)

install(TARGETS ttkManifoldCheckGui RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})
