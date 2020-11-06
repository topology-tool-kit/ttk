set(TTK_INSTALL_BINARY_DIR "bin")

set(cmd_dirs "*/cmd")
set(gui_dirs "*/gui")

# still build cmd if no rendering
if(NOT TARGET ttkUserInterfaceBase)
  message(WARNING "ttkUserInterfaceBase is disabled, no gui standalone")
  set(gui_dirs "")
endif()

# warning: use deprecated `add_compile_option` here to avoid modifying
# every standalone CMakeLists.txt
add_compile_options(${TTK_COMPILER_FLAGS})

file(GLOB STANDALONE_DIRS ${cmd_dirs} ${gui_dirs} "./*")
foreach(STANDALONE ${STANDALONE_DIRS})
  if(IS_DIRECTORY ${STANDALONE} AND EXISTS ${STANDALONE}/CMakeLists.txt)
    add_subdirectory(${STANDALONE})
  endif()
endforeach()

install(DIRECTORY textures DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/ttk")
