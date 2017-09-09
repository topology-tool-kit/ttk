cmake_minimum_required(VERSION 3.2)

project(morseSmaleComplexCmd)
add_executable(morseSmaleComplexCmd main.cpp)
target_link_libraries(morseSmaleComplexCmd ttkMorseSmaleComplex ttkProgramBase)
install(TARGETS morseSmaleComplexCmd RUNTIME DESTINATION bin/ttk)

