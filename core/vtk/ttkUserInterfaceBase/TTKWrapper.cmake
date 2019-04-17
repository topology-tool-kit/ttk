ttk_register_vtk_filter()

# TODO test
target_compile_definitions(ttkUserInterfaceBase
  PRIVATE
  TTK_INSTALL_ASSETS_DIR="${CMAKE_INSTALL_PREFIX}/share/ttk"
  )
