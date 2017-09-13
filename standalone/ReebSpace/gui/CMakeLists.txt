cmake_minimum_required(VERSION 3.2)

project(reebSpaceGui)
add_executable(reebSpaceGui main.cpp)
target_link_libraries(reebSpaceGui ttkReebSpace ttkUserInterfaceBase)

install(TARGETS reebSpaceGui RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

