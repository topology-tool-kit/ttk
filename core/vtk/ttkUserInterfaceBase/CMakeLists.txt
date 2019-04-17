ttk_add_vtk_library(ttkUserInterfaceBase
  SOURCES
    ttkUserInterfaceBase.cpp
  HEADERS
    ttkUserInterfaceBase.h
  LINK
    ttkProgramBase
    ttkTextureMapFromField
    )

target_compile_definitions(ttkUserInterfaceBase
  PRIVATE
    TTK_INSTALL_ASSETS_DIR="${CMAKE_INSTALL_PREFIX}/share/ttk"
    )
