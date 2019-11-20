ttk_add_vtk_module()

if(NOT TTK_PYTHON_MAJOR_VERSION)
  set(VTK_MODULE_ENABLE_ttkDimensionReduction "NO" CACHE BOOL "Disable dimension reduction" FORCE)
  message(STATUS "No python found, disable ttk dimension reduction filter")
endif()
