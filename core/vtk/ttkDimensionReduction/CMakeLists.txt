ttk_add_vtk_module()

if(NOT PYTHON_INCLUDE_DIRS)
  set(VTK_MODULE_ENABLE_ttkDimensionReduction "NO" CACHE BOOL "Disable dimension reduction" FORCE)
  message(STATUS "No python found, disable ttk dimension reduction filter")
endif()
