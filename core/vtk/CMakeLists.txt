include(functions.cmake)

option(TTK_FORCE_VTK_PARAVIEW "Force to use the VTK bundled from ParaView" ON)
mark_as_advanced(TTK_FORCE_VTK_PARAVIEW)

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

include_directories(${VTK_INCLUDE_DIRS})

set(VTKWRAPPER_LIB_LIST "")
file(GLOB VTKWRAPPER_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ttk*)
foreach(VTKWRAPPER ${VTKWRAPPER_DIRS})
  if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${VTKWRAPPER})
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${VTKWRAPPER})
    set(VTKWRAPPER_LIB_LIST "${VTKWRAPPER};${VTKWRAPPER_LIB_LIST}")
  endif()
endforeach()

# Add a bundle library that will just link all the TTK vtkWrappers
add_library(ttkAll INTERFACE)
target_include_directories(ttkAll INTERFACE ${VTK_INCLUDE_DIRS})
target_link_libraries(ttkAll INTERFACE ${VTKWRAPPER_LIB_LIST})

install(TARGETS ttkAll
  EXPORT TTKVTKTargets
  RUNTIME DESTINATION bin/ttk
  ARCHIVE DESTINATION lib/ttk
  LIBRARY DESTINATION lib/ttk
  )

install(EXPORT TTKVTKTargets
  NAMESPACE ttk::vtk::
  DESTINATION lib/cmake/ttk
  )

configure_file(TTKVTKConfig.cmake.in TTKVTKConfig.cmake @ONLY)
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/TTKVTKConfig.cmake"
  DESTINATION lib/cmake/ttk
  )

set(_build_vtk_python_default OFF)
if (VTK_WRAP_PYTHON AND (${VTK_MAJOR_VERSION} EQUAL 8) AND (VTK_VERSION VERSION_LESS 8.90.0))
  set(_build_vtk_python_default ON)
endif()
option(TTK_BUILD_VTK_PYTHON_MODULE "Build VTK Python Module" ${_build_vtk_python_default})

if (TTK_BUILD_VTK_PYTHON_MODULE)
  if (NOT ((${VTK_MAJOR_VERSION} EQUAL 8) AND (VTK_VERSION VERSION_LESS 8.90.0)))
    message(FATAL_ERROR "Building VTK Python module is only supported for VTK version 8.")
  endif()
  if (NOT VTK_WRAP_PYTHON)
    message(FATAL_ERROR "VTK needs to be built with Python wrapping enabled in order to build the VTK Python module.")
  endif()

  # name of the python module
  set(module topologytoolkit)
  # VTK dependencies
  set(module_depends 
    vtkCommonExecutionModel
    vtkIOParallelXML)

  # find VTK prefix
  if (VTK_INSTALL_PREFIX)
    set(_vtk_path ${VTK_INSTALL_PREFIX})
  elseif (VTK_DIR)
    set(_vtk_path ${VTK_DIR})
  else ()
    set(_vtk_path ${ParaView_DIR})
  endif()
  if (_vtk_path AND VTK_PYTHONPATH)
    # guess relative python module path from VTK config
    if (NOT IS_ABSOLUTE ${VTK_PYTHONPATH})
      file(RELATIVE_PATH _python_module_dir ${_vtk_path} ${_vtk_path}/${VTK_PYTHONPATH})
    else ()
      file(RELATIVE_PATH _python_module_dir ${_vtk_path} ${VTK_PYTHONPATH})
    endif ()
  endif()
  set(TTK_PYTHON_MODULE_DIR
    "${_python_module_dir}"
    CACHE STRING
    "Subdirectory where VTK Python module will be installed")

  # get headers and include directories of vtkwrappers
  set(headers)
  set(include_dirs)
  foreach (vtkwrapper ${VTKWRAPPER_LIB_LIST})
    get_target_property(_source_dir ${vtkwrapper} SOURCE_DIR)
    get_target_property(_headers_relative ${vtkwrapper} TTK_HEADERS)
    set(include_dirs ${include_dirs} ${_source_dir})
    foreach (header ${_headers_relative})
      if (IS_ABSOLUTE ${header})
        set(headers ${headers} "${header}")
      else ()
        set(headers ${headers} "${_source_dir}/${header}")
      endif ()  
    endforeach ()
  endforeach ()

  # add dummy library (needed for vtk_wrap_hierarchy)
  add_library(${module} EXCLUDE_FROM_ALL ${headers})
  set_target_properties(${module} PROPERTIES LINKER_LANGUAGE CXX)
  target_include_directories(${module} PUBLIC ${include_dirs})

  # generate ttkHierarchy.txt
  # only those classes, that have a full inheritance chain that ends in vtkObject will be wrapped
  set(MODULE_HIERARCHY_NAME ${module}Hierarchy)
  set(${module}_WRAP_DEPENDS ${module_depends})
  include(vtkWrapHierarchy)
  vtk_wrap_hierarchy(${module} ${CMAKE_CURRENT_BINARY_DIR} "${headers}")
  set(KIT_HIERARCHY_FILE ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.txt)
  set(LIB_HIERARCHY_STAMP ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.stamp.txt)

  add_library(${module}Python MODULE ${module}PythonInit.cxx)
  set_target_properties(${module}Python PROPERTIES PREFIX "")
  if (WIN32 AND NOT CYGWIN)
    set_target_properties(${module}Python PROPERTIES SUFFIX ".pyd")
  endif ()  
  target_include_directories(${module}Python PUBLIC ${include_dirs})

  include(vtkWrapPython)
  vtk_wrap_python3(${module}Python ${module}Python_SOURCES "${headers}")

  add_library(${module}PythonD STATIC ${${module}Python_SOURCES})
  set_target_properties(${module}PythonD PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  target_include_directories(${module}PythonD PUBLIC ${include_dirs})
  target_link_libraries(${module}PythonD ${VTKWRAPPER_LIB_LIST})
  target_link_libraries(${module}PythonD
    ${VTK_LIBRARIES}
    vtkWrappingPythonCore
    ${VTK_PYTHON_LIBRARIES})

  # link against VTK PythonD modules
  set(_pythond_libs)
  foreach(vtk_module ${module_depends})
    set(_pythond_libs ${_pythond_libs} ${vtk_module}PythonD)
  endforeach()
  target_link_libraries(${module}Python ${module}PythonD ${_pythond_libs})

  set_target_properties(${module}Python PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TTK_PYTHON_MODULE_DIR}/${module})
  install(TARGETS ${module}Python DESTINATION ${TTK_PYTHON_MODULE_DIR}/${module})
  # set rpath such that vtkwrappers can be found
  set_target_properties(${module}Python PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")

  # generate __init__.py
  set(module_py "from __future__ import absolute_import\nfrom .${module}Python import *\n")
  file(GENERATE
    OUTPUT "${CMAKE_BINARY_DIR}/${TTK_PYTHON_MODULE_DIR}/${module}/__init__.py"
    CONTENT "${module_py}")
  install(
    FILES       "${CMAKE_BINARY_DIR}/${TTK_PYTHON_MODULE_DIR}/${module}/__init__.py"
    DESTINATION "${TTK_PYTHON_MODULE_DIR}/${module}")
endif()
