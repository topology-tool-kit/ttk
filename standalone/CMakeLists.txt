set(TTK_INSTALL_BINARY_DIR "bin")

include(../core/vtk/ttkProgramBase/TTKWrapper.cmake)
include(../core/vtk/ttkUserInterfaceBase/TTKWrapper.cmake)

file(GLOB STANDALONE_DIRS */cmd */gui)
foreach(STANDALONE ${STANDALONE_DIRS})
  if(IS_DIRECTORY ${STANDALONE})
    add_subdirectory(${STANDALONE})
  endif()
endforeach()

install(DIRECTORY textures DESTINATION share/ttk)
