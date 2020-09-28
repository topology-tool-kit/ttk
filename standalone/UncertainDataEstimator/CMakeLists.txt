cmake_minimum_required(VERSION 3.2)

project(ttkUncertainDataEstimatorCmd)

if(TARGET ttkUncertainDataEstimator)
  add_executable(${PROJECT_NAME} main.cpp Editor.cpp)
  target_link_libraries(${PROJECT_NAME}
    PRIVATE
      ttkUncertainDataEstimator
      VTK::IOImage
      VTK::IOLegacy
      VTK::IOXML
    )
  set_target_properties(${PROJECT_NAME}
    PROPERTIES
      INSTALL_RPATH
        "${CMAKE_INSTALL_RPATH}"
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
