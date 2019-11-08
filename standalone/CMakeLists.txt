set(TTK_INSTALL_BINARY_DIR "bin")

set(cmd_dirs "*/cmd")
set(gui_dirs "*/gui")

# still build cmd if no rendering
if(NOT TARGET ttkUserInterfaceBase)
  message(WARNING "ttkUserInterfaceBase is disabled, no gui standalone")
  set(gui_dirs "")
endif()

# TODO investigate (mac)
# set(CMAKE_SKIP_BUILD_RPATH TRUE)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

file(GLOB STANDALONE_DIRS ${cmd_dirs} ${gui_dirs})
foreach(STANDALONE ${STANDALONE_DIRS})
  if(IS_DIRECTORY ${STANDALONE})
    add_subdirectory(${STANDALONE})
  endif()
endforeach()

install(DIRECTORY textures DESTINATION share/ttk)
