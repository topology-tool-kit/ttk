include(functions.cmake)

# build VTK wrappers explicitely and export them
file(GLOB VTKWRAPPER_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ttk*)
foreach(VTKWRAPPER ${VTKWRAPPER_DIRS})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${VTKWRAPPER}/TTKWrapper.cmake)
    include(${VTKWRAPPER}/TTKWrapper.cmake)
    set(TTK_INSTALL_TARGETS "${VTKWRAPPER};${TTK_INSTALL_TARGETS}")
  endif()
endforeach()

# Add a bundle library that will just link all the TTK vtkWrappers
add_library(ttkAll INTERFACE)
target_include_directories(ttkAll INTERFACE ${VTK_INCLUDE_DIRS})
target_link_libraries(ttkAll INTERFACE ${TTK_INSTALL_TARGETS})

message(STATUS "TTK modules: ${TTK_INSTALL_TARGETS}")

vtk_module_autoinit(
  TARGETS
    ttkAll
  MODULES
    ${VTK_LIBRARIES}
  )

install(
  TARGETS
    ttkAll
  EXPORT
    TTKVTKTargets
  RUNTIME DESTINATION
    bin/ttk
  ARCHIVE DESTINATION
    lib/ttk
  LIBRARY DESTINATION
    lib/ttk
  )

install(
  EXPORT
    TTKVTKTargets
  NAMESPACE
    ttk::vtk::
  DESTINATION
    lib/cmake/ttk
    )

configure_file(TTKVTKConfig.cmake.in TTKVTKConfig.cmake @ONLY)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/TTKVTKConfig.cmake"
  DESTINATION
    lib/cmake/ttk
    )

