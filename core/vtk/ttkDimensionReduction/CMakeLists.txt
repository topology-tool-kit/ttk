ttk_add_vtk_module()

if(NOT TTK_ENABLE_SCIKIT_LEARN)
  set(VTK_MODULE_ENABLE_ttkDimensionReduction "NO"
    CACHE STRING "Enable the ttkDimensionReduction module." FORCE)
  message(STATUS "No Python or scikit-learn found, ttkDimensionReduction filter disabled")
endif()
