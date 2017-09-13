cmake_minimum_required(VERSION 3.2)

project(morseSmaleComplexGui)
add_executable(morseSmaleComplexGui main.cpp)
target_link_libraries(morseSmaleComplexGui ttkMorseSmaleComplex ttkUserInterfaceBase)

install(TARGETS morseSmaleComplexGui RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

