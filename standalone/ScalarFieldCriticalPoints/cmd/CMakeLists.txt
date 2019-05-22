cmake_minimum_required(VERSION 3.2)

project(ttkScalarFieldCriticalPointsCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(NOT DEFINED TTK_BUILD_SCALARFIELDCRITICALPOINTS_FILTER)
  # Allows to disable each filter
  option(TTK_BUILD_SCALARFIELDCRITICALPOINTS_FILTER "Build the filter ScalarFieldCriticalPoints" ON)
  mark_as_advanced(TTK_BUILD_SCALARFIELDCRITICALPOINTS_FILTER)
endif()

if(${TTK_BUILD_SCALARFIELDCRITICALPOINTS_FILTER})
  add_executable(ttkScalarFieldCriticalPointsCmd main.cpp)
  target_link_libraries(ttkScalarFieldCriticalPointsCmd
    PRIVATE
      ttkScalarFieldCriticalPoints
      ttkProgramBase
      )
  install(
    TARGETS
      ${PROJECT_NAME}
    RUNTIME DESTINATION
      ${TTK_INSTALL_BINARY_DIR}
      )
endif()

