cmake_minimum_required(VERSION 3.2)

project(blankGui)

add_executable(blankGui main.cpp)
target_link_libraries(blankGui ${VTK_LIBRARIES}
	ttkBlank ttkProgramBase ttkUserInterfaceBase)

install(TARGETS blankGui RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})
install(DIRECTORY textures DESTINATION ${TTK_INSTALL_BINARY_DIR})

