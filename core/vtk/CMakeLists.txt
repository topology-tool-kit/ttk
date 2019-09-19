vtk_module_find_modules(
  TTK_VTK_MODULE_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}/${VTKWRAPPER}"
  )

vtk_module_scan(
  MODULE_FILES
    ${TTK_VTK_MODULE_FILE}
  PROVIDES_MODULES
    TTK_ENABLED_MODULES
  WANT_BY_DEFAULT
    ON
  )

if(TTK_BUILD_PARAVIEW_PLUGINS)
  file(GLOB VTKWRAPPER_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ttk*)
  foreach(VTKWRAPPER ${VTKWRAPPER_DIRS})
    ttk_register_vtk_filter(${CMAKE_CURRENT_SOURCE_DIR}/${VTKWRAPPER})
  endforeach()
else()
  vtk_module_build(
    MODULES
      ${TTK_ENABLED_MODULES}
    )
endif()

# Add a bundle library that will just link all the TTK vtkWrappers
add_library(ttkAll INTERFACE)

vtk_module_autoinit(
  TARGETS
    ttkAll
  MODULES
    ${TTK_ENABLED_MODULES}
  )

if(NOT TTK_BUILD_PARAVIEW_PLUGINS)
  # when building against VTK, we need the modules
  # in addition to ttkAll.
  set(TTK_EXPORT_TARGETS ${TTK_BUILD_PARAVIEW_PLUGINS})
endif()

# Use folder as configured by the distribution
include(GNUInstallDirs)
install(
  TARGETS
    ttkAll
    ${TTK_EXPORT_TARGETS}
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
