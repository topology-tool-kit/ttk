cmake_minimum_required(VERSION 3.2)

project(ttkBlankGui)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(NOT DEFINED TTK_BUILD_BLANK_FILTER)
  # Allows to disable each filter
  option(TTK_BUILD_BLANK_FILTER "Build the filter Blank" ON)
  mark_as_advanced(TTK_BUILD_BLANK_FILTER)
endif()

if(${TTK_BUILD_BLANK_FILTER})
  add_executable(ttkBlankGui main.cpp)
  target_link_libraries(ttkBlankGui
    PRIVATE
      ttkBlank
      ttkProgramBase
      ttkUserInterfaceBase
      )
  install(
    TARGETS
      ${PROJECT_NAME}
    RUNTIME DESTINATION
      ${TTK_INSTALL_BINARY_DIR}
      )
endif()
