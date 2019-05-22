cmake_minimum_required(VERSION 3.2)

project(ttkManifoldCheckCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(NOT DEFINED TTK_BUILD_MANIFOLDCHECK_FILTER)
  # Allows to disable each filter
  option(TTK_BUILD_MANIFOLDCHECK_FILTER "Build the filter ManifoldCheck" ON)
  mark_as_advanced(TTK_BUILD_MANIFOLDCHECK_FILTER)
endif()

if(${TTK_BUILD_MANIFOLDCHECK_FILTER})
  add_executable(ttkManifoldCheckCmd main.cpp)
  target_link_libraries(ttkManifoldCheckCmd
    PRIVATE
      ttkManifoldCheck
      ttkProgramBase
      )
  install(
    TARGETS
      ${PROJECT_NAME}
    RUNTIME DESTINATION
      ${TTK_INSTALL_BINARY_DIR}
      )
endif()

