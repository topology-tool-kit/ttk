ttk_add_vtk_module()

if(NOT TTK_ENABLE_SCIKIT_LEARN)
  set(VTK_MODULE_ENABLE_ttkDimensionReduction "NO" CACHE BOOL "Disable dimension reduction" FORCE)
  message(STATUS "No python or scikit-learn found, disable ttk dimension reduction filter")
endif()
