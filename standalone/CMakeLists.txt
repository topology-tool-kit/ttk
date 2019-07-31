# find VTK, if none provided but paraview is given,
# use the VTK provided by paraview with a message
if(NOT ${TTK_FORCE_VTK_PARAVIEW})
  find_package(VTK QUIET)
endif()
if(NOT DEFINED VTK_USE_FILE)
  find_package(ParaView)
  if(NOT DEFINED VTK_USE_FILE)
    # VTK error message
    find_package(VTK REQUIRED)
  else()
    message(STATUS "No VTK provided, fallback to ParaView")
  endif()
endif()

if(DEFINED VTK_USE_FILE)
    include(${VTK_USE_FILE})
endif()

set(TTK_INSTALL_BINARY_DIR "bin")
file(GLOB STANDALONE_DIRS */cmd */gui)
foreach(STANDALONE ${STANDALONE_DIRS})
  if (IS_DIRECTORY ${STANDALONE})
    add_subdirectory(${STANDALONE})
  endif()
endforeach()

install(DIRECTORY textures DESTINATION share/ttk)
