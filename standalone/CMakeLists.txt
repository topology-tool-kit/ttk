set(TTK_INSTALL_BINARY_DIR "bin")

if(NOT TARGET ttkProgramBase)
  include(../core/vtk/ttkProgramBase/TTKWrapper.cmake)
endif()
if(NOT TARGET ttkUserInterfaceBase)
  include(../core/vtk/ttkUserInterfaceBase/TTKWrapper.cmake)
endif()

file(GLOB STANDALONE_DIRS */cmd */gui)
foreach(STANDALONE ${STANDALONE_DIRS})
  if(IS_DIRECTORY ${STANDALONE})
    add_subdirectory(${STANDALONE})
  endif()
endforeach()

install(DIRECTORY textures DESTINATION share/ttk)
