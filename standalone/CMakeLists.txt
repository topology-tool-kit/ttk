set(TTK_INSTALL_BINARY_DIR "bin")

if(NOT TARGET ttkProgramBase)
  ttk_register_vtk_filter("../core/vtk/ttkProgramBase")
endif()
if(NOT TARGET ttkUserInterfaceBase)
  ttk_register_vtk_filter("../core/vtk/ttkUserInterfaceBase")
endif()

# TODO investigate (mac)
# set(CMAKE_SKIP_BUILD_RPATH TRUE)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

file(GLOB STANDALONE_DIRS */cmd */gui)
foreach(STANDALONE ${STANDALONE_DIRS})
  if(IS_DIRECTORY ${STANDALONE})
    add_subdirectory(${STANDALONE})
  endif()
endforeach()

install(DIRECTORY textures DESTINATION share/ttk)
