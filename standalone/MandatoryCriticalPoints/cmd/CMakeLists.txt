cmake_minimum_required(VERSION 3.2)

project(ttkMandatoryCriticalPointsCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(NOT DEFINED TTK_BUILD_MANDATORYCRITICALPOINTS_FILTER)
  # Allows to disable each filter
  option(TTK_BUILD_MANDATORYCRITICALPOINTS_FILTER "Build the filter MandatoryCriticalPoints" ON)
  mark_as_advanced(TTK_BUILD_MANDATORYCRITICALPOINTS_FILTER)
endif()

if(${TTK_BUILD_MANDATORYCRITICALPOINTS_FILTER})
  add_executable(ttkMandatoryCriticalPointsCmd main.cpp)
  target_link_libraries(ttkMandatoryCriticalPointsCmd
    PRIVATE
      ttkMandatoryCriticalPoints
      ttkProgramBase
      )
  install(
    TARGETS
      ${PROJECT_NAME}
    RUNTIME DESTINATION
      ${TTK_INSTALL_BINARY_DIR}
      )
endif()

