cmake_minimum_required(VERSION 3.2)

project(ttkUncertainDataEstimatorCmd)

if(NOT DEFINED TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER)
  # Allows to disable each filter
  option(TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER "Build the filter UncertainDataEstimator" ON)
  mark_as_advanced(TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER)
endif()

if(${TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER})
  add_executable(${PROJECT_NAME} main.cpp Editor.cpp)
  target_link_libraries(${PROJECT_NAME}
    PRIVATE
      ttkUncertainDataEstimator
      ttkProgramBase
      )
  set_target_properties(${PROJECT_NAME}
    PROPERTIES
      INSTALL_RPATH
        "${CMAKE_INSTALL_PREFIX}/lib64/TopologyToolKit;${CMAKE_INSTALL_PREFIX}/lib/ttk")
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

