vtk_module_find_modules(
  TTK_VTK_MODULE_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}/${VTKWRAPPER}"
  )

vtk_module_scan(
  MODULE_FILES
    ${TTK_VTK_MODULE_FILE}
  PROVIDES_MODULES
    TTK_PROVIDED_MODULES
  WANT_BY_DEFAULT
    ON
  )

set(TTK_ENABLED_MODULES "")

# Whitelist mechanism
foreach(TTK_MODULE ${TTK_PROVIDED_MODULES})
  ttk_get_target(${TTK_MODULE} TTK_TARGET)
  if(NOT "${TTK_TARGET}" STREQUAL "")
    list(APPEND TTK_ENABLED_MODULES ${TTK_TARGET})
  endif()
endforeach()

if(NOT ${TTK_BUILD_PARAVIEW_PLUGINS})
  # when build against ParaView,
  # vtk_module_build is already called
  # by the paraview_add_plugin
  vtk_module_build(
    MODULES
      ${TTK_ENABLED_MODULES}
    )
endif()

if(VTK_WRAP_PYTHON)
  set(TTK_BUILD_VTK_PYTHON_MODULE
    ON
    CACHE
    BOOL
    "Build the python module"
    )
  if(NOT DEFINED TTK_PYTHON_MODULE_DIR)
    vtk_module_python_default_destination(PYTHON_SITE_PACKAGES_SUFFIX)
    set(TTK_PYTHON_MODULE_DIR
      ${PYTHON_SITE_PACKAGES_SUFFIX}
      CACHE
      PATH
      "Where the TTK python module in installed"
      )
  endif()
  vtk_module_wrap_python(
    MODULES
      ${TTK_ENABLED_MODULES}
    WRAPPED_MODULES
      TTK_WRAPPED_MODULES
    INSTALL_EXPORT
      TTKPython
    PYTHON_PACKAGE
      "topologytoolkit"
    MODULE_DESTINATION
      "${TTK_PYTHON_MODULE_DIR}"
    CMAKE_DESTINATION
      "${CMAKE_INSTALL_LIBDIR}/ttk/python"
    LIBRARY_DESTINATION
      "${CMAKE_INSTALL_LIBDIR}"
    INSTALL_HEADERS
      ON
    TARGET
      TTK::Python
    )
  export(
    EXPORT
      TTKPython
    NAMESPACE
      TTK::
    FILE
      "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/ttk/python/TTKPython-targets.cmake"
    )
  install(
    EXPORT
      TTKPython
    NAMESPACE
      TTK::
    FILE
      TTKPython-targets.cmake
    DESTINATION
      "${CMAKE_INSTALL_LIBDIR}/ttk/python"
    )
else()
  set(TTK_BUILD_VTK_PYTHON_MODULE
    OFF
    CACHE
    BOOL
    "Build the python module"
    )
endif()

# Add a bundle library that will just link all the TTK vtkWrappers
add_library(ttkTMP INTERFACE)

vtk_module_autoinit(
  TARGETS
    ttkTMP
  MODULES
    ${TTK_ENABLED_MODULES}
  )

set(TTK_ENABLED_MODULES ${TTK_ENABLED_MODULES} CACHE INTERNAL "")
