vtk_module_find_modules(
  TTK_MODULE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/${VTKWRAPPER}"
  )

vtk_module_scan(
  MODULE_FILES
    ${TTK_MODULE_FILES}
  PROVIDES_MODULES
    TTK_ENABLED_MODULES
  WANT_BY_DEFAULT
    ON
  )

if (TTK_BUILD_PARAVIEW_PLUGINS)
  file(GLOB VTKWRAPPER_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ttk*)
  foreach(VTKWRAPPER ${VTKWRAPPER_DIRS})
    ttk_register_vtk_filter(${CMAKE_CURRENT_SOURCE_DIR}/${VTKWRAPPER})
    list(APPEND TTK_ENABLED_MODULES ${VTKWRAPPER})
  endforeach()
else()
  vtk_module_build(
    MODULES
      ${TTK_ENABLED_MODULES}
    )
endif()

# Add a bundle library that will just link all the TTK vtkWrappers
add_library(ttkAll INTERFACE)
target_include_directories(ttkAll INTERFACE ${VTK_INCLUDE_DIRS})

vtk_module_autoinit(
  TARGETS
    ttkAll
  MODULES
    ${TTK_ENABLED_MODULES}
  )

# Use folder as configured by the distribution
include(GNUInstallDirs)
install(
  TARGETS
    ttkAll
  EXPORT
    TTKVTKTargets
  RUNTIME DESTINATION
    ${CMAKE_INSTALL_BINDIR}/ttk
  ARCHIVE DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/ttk
  LIBRARY DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/ttk
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
