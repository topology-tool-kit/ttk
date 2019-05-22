cmake_minimum_required(VERSION 3.2)

project(ttkUncertainDataEstimatorCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(NOT DEFINED TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER)
  # Allows to disable each filter
  option(TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER "Build the filter UncertainDataEstimator" ON)
  mark_as_advanced(TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER)
endif()

if(${TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER})
  add_executable(ttkUncertainDataEstimatorCmd main.cpp Editor.cpp)
  target_link_libraries(ttkUncertainDataEstimatorCmd
    PRIVATE
      ttkUncertainDataEstimator
      ttkProgramBase
      )
  target_include_directories(ttkUncertainDataEstimatorCmd
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )
  install(
    TARGETS
      ${PROJECT_NAME}
    RUNTIME DESTINATION
      ${TTK_INSTALL_BINARY_DIR}
      )
endif()

