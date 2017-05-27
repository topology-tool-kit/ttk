# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

if(NOT TTK_PKG)

  # init project global variables
  set(PROJECT_DEP "" 
    CACHE INTERNAL "PROJECT_DEP")
  set(PROJECT_FLAGS "" 
    CACHE INTERNAL "PROJECT_FLAGS")
  set(PROJECT_SRC "" 
    CACHE INTERNAL "PROJECT_SRC")

  foreach(LIB ${LIB_LIST})
    set(${LIB} "" CACHE INTERNAL "${LIB}")
  endforeach(LIB)
  set(LIB_LIST "" CACHE INTERNAL "LIB_LIST")

  # specify the path to the core packages
  find_path(BASECODE_DIR baseCode.cmake
    PATHS
      ${TTK_DIR}/baseCode/
      baseCode/
      ../baseCode/
      ../../baseCode/
  )
  include(${BASECODE_DIR}/baseCode.cmake)

  # cmake helper functions
  function(ttk_add_baseCode_package package)
    if(NOT ${package}_PATH)
      ttk_find_package(${package})
      if(NOT ${package}_PATH)
        message(FATAL_ERROR 
          "Could not find the package '${package}.cmake'!") 
      endif(NOT ${package}_PATH)
    endif(NOT ${package}_PATH) 
  endfunction(ttk_add_baseCode_package)

  function(ttk_add_cflags flags)
    set(PROJECT_FLAGS ${PROJECT_FLAGS}
      ${flags}
      CACHE INTERNAL "PROJECT_FLAGS")
  endfunction(ttk_add_cflags)

  function(ttk_add_dep dep)
    set(PROJECT_DEP 
      ${PROJECT_DEP} 
      ${dep} 
      CACHE INTERNAL "PROJECT_DEP")
  endfunction(ttk_add_dep)

  function(ttk_add_external_package header lib)
    option(with${lib} "Enable ${lib} support" true)

    if(with${lib})
      find_path(${header}_PATH ${header}
        PATHS
          ${ARGN}
      )
      if(NOT ${header}_PATH)
        message(STATUS
"ttk -------------------------------------------------------------------------")
        message(STATUS "Header ${header} not found. External package disabled.")
        message(STATUS
"-----------------------------------------------------------------------------")
        set(with${lib} OFF)
      endif(NOT ${header}_PATH)

      if(${header}_PATH)
        include_directories(${${header}_PATH})
        find_library(${lib}_PATH 
          NAMES ${lib}
          PATHS
            ${ARGN}
        )
        if(NOT ${lib}_PATH)
          message(STATUS
"ttk -------------------------------------------------------------------------")
          message(STATUS "Library ${lib} not found. External package disabled.")
          message(STATUS
"-----------------------------------------------------------------------------")
          set(with${lib} OFF)
        endif(NOT ${lib}_PATH)

        if(${lib}_PATH)
          ttk_add_option(with${lib})
          ttk_add_dep(${${lib}_PATH})
        endif(${lib}_PATH)
      endif(${header}_PATH)
    endif(with${lib})

    message(STATUS
"ttk -------------------------------------------------------------------------")
    message(STATUS
    "${lib} support: ${with${lib}}   (-Dwith${lib}=)")
    message(STATUS
"-----------------------------------------------------------------------------")
  endfunction(ttk_add_external_package)

  function(ttk_add_option option)
    set(PROJECT_FLAGS "${PROJECT_FLAGS} -D${option}"
      CACHE INTERNAL "PROJECT_FLAGS")
  endfunction(ttk_add_option)

  function(ttk_add_optional_baseCode_package package)
    option(with${package} "Enable ${package} support" true)

    if(${with${package}})
      if(NOT ${package}_PATH)
        ttk_find_package(${package})
        if(${package}_PATH)
          ttk_add_option(with${package})
        elseif(NOT ${package}_PATH)
          set(with${package} OFF)
          message(STATUS
"ttk -------------------------------------------------------------------------")
          message(STATUS
            "Package ${package} not found. Option disabled.")
          message(STATUS
"-----------------------------------------------------------------------------")
        endif(${package}_PATH)
      endif(NOT ${package}_PATH) 
    endif(${with${package}})

    message(STATUS
"ttk -------------------------------------------------------------------------")
    message(STATUS
  "${package} support: ${with${package}}   (-Dwith${package}=)")
    message(STATUS
"-----------------------------------------------------------------------------")
  endfunction(ttk_add_optional_baseCode_package)

  function(ttk_add_optional_vtkWrapper_package package)
    if(NOT withVTK)
      message(FATAL_ERROR "Trying to compile a VTK wrapper with VTK disabled!")
    endif(NOT withVTK)
    ttk_add_optional_baseCode_package(${package})
  endfunction(ttk_add_optional_vtkWrapper_package)

  function(ttk_add_source source)
    include_directories(${CMAKE_CURRENT_LIST_DIR})
    set(PROJECT_SRC ${PROJECT_SRC} 
      ${CMAKE_CURRENT_LIST_DIR}/${source}
      CACHE INTERNAL "PROJECT_SRC")
  endfunction(ttk_add_source)

  function(ttk_add_vtkWrapper_package package)
    if(NOT withVTK)
      message(FATAL_ERROR "Trying to compile a VTK wrapper with VTK disabled!")
    endif(NOT withVTK)
    ttk_add_baseCode_package(${package})
  endfunction(ttk_add_vtkWrapper_package)

  function(ttk_find_package package)
    set(${package}_PATH ${package}_PATH-NOTFOUND 
      CACHE INTERNAL "${package}_PATH")
    find_path(${package}_PATH "${package}.cmake"
      PATHS
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../../sandbox/baseCode/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../../sandbox/vtkWrappers/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../sandbox/baseCode/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../sandbox/vtkWrappers/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../../sandbox/baseCode/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../../sandbox/vtkWrappers/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../sandbox/baseCode/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../sandbox/vtkWrappers/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/sandbox/baseCode/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/sandbox/vtkWrappers/${package}/
        ${TTK_DIR}/baseCode/${package}/
        ${TTK_DIR}/vtkWrappers/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../../sandbox/*/sandbox/baseCode/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../sandbox/*/sandbox/baseCode/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../../sandbox/*/sandbox/baseCode/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/../sandbox/*/sandbox/baseCode/${package}/
        ${CMAKE_CURRENT_SOURCE_DIR}/sandbox/*/sandbox/baseCode/${package}/

      NO_DEFAULT_PATH)
    if(${package}_PATH)
      # the package may just be a template (no call to add_source)
      include_directories("${${package}_PATH}")
      include("${${package}_PATH}/${package}.cmake")
      set(PACKAGE_LIST ${PACKAGE_LIST} ${package}
        CACHE INTERNAL "PACKAGE_LIST")
    endif(${package}_PATH)
  endfunction(ttk_find_package)

  function(ttk_wrapup_binary project)
    if(APPLE)
      add_executable(${project} MACOSX_BUNDLE ${PROJECT_SRC})
    else()
      add_executable(${project} ${PROJECT_SRC})
    endif(APPLE)
    ttk_wrapup_flags(${project})
    install(
      TARGETS ${project}
      DESTINATION ${TTK_BINARY_INSTALL_DIR})
  endfunction(ttk_wrapup_binary)

  function(ttk_wrapup_flags project)
    
    # set up compilation flags pulled out by the called modules 
    set_target_properties(${project}
      PROPERTIES
      COMPILE_FLAGS 
      "${PROJECT_FLAGS}")
  
    # specify the required libraries for linkage, as pulled out by the called 
    # modules (${PROJECT_DEP})
    target_link_libraries(${project} "${PROJECT_DEP}")

    message(STATUS
"ttk -------------------------------------------------------------------------")
    message(STATUS "ttk project: ${project}")
    get_property(PROJECT_INCLUDES DIRECTORY 
      ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
    message(STATUS "ttk includes: ${PROJECT_INCLUDES}")
    message(STATUS "ttk deps: ${PROJECT_DEP}")
    message(STATUS "ttk flags: ${PROJECT_FLAGS}")
    message(STATUS "ttk src: ${PROJECT_SRC}")
    message(STATUS "ttk packages: ${PACKAGE_LIST}")
    message(STATUS
"-----------------------------------------------------------------------------")
  # clear the package markers
    foreach(package ${PACKAGE_LIST})
      set(${package}_PATH ${package}_PATH-NOTFOUND 
        CACHE INTERNAL "${package}_PATH")
    endforeach(package)
    set(PACKAGE_LIST ""
      CACHE INTERNAL "PACKAGE_LIST")
  endfunction(ttk_wrapup_flags)

  function(ttk_wrapup_library library source)

    # link only once!
    foreach(LIB ${LIB_LIST})
      if("${LIB}" MATCHES " ${source}")
        return()
      endif("${LIB}" MATCHES " ${source}")
    endforeach(LIB)

    string(REPLACE " " ";" sourceList ${source})

    foreach(sourceFile ${sourceList})
      if("${TTK_BUILD_MODE}" MATCHES "ParaView")
        add_library(${library} OBJECT 
          ${CMAKE_CURRENT_LIST_DIR}/${sourceFile})
      elseif(NOT "${TTK_BUILD_MODE}" MATCHES "ParaView")
        ttk_add_source(${sourceFile})
      endif("${TTK_BUILD_MODE}" MATCHES "ParaView")
    endforeach(sourceFile)

    if("${TTK_BUILD_MODE}" MATCHES "ParaView")
      # use the project flags
      set_target_properties(${library}
        PROPERTIES
        COMPILE_FLAGS
        "${PROJECT_FLAGS}")
      add_library(${PROJECT_NAME} SHARED $<TARGET_OBJECTS:${library}>)
    endif("${TTK_BUILD_MODE}" MATCHES "ParaView")

    set(LIB_LIST "${LIB_LIST} ${source}" CACHE INTERNAL "LIB_LIST")

  endfunction(ttk_wrapup_library)

  function(ttk_wrapup_paraview_plugin
    plugin_name
    plugin_version)

    if("${TTK_BUILD_MODE}" MATCHES "ParaView")
      
      add_paraview_plugin(
        # name
        ${plugin_name}
        # version
        "${plugin_version}"                         
        SERVER_MANAGER_XML "${plugin_name}.xml"
        SERVER_MANAGER_SOURCES "${PROJECT_SRC}"
      )
    endif("${TTK_BUILD_MODE}" MATCHES "ParaView")

    ttk_wrapup_flags(${plugin_name})
    install(
      TARGETS ${plugin_name}
      DESTINATION ${TTK_PLUGIN_INSTALL_DIR})
  endfunction(ttk_wrapup_paraview_plugin)

  # message about the location of the common code base
  message(STATUS
"ttk -------------------------------------------------------------------------")
  message(STATUS "ttk path: ${TTK_DIR}")
  message(STATUS
"-----------------------------------------------------------------------------")

  set(withVTK false)

  set(TTK_BINARY_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin"
    CACHE PATH "Directory where TTK binary programs will be installed")
  set(CMAKE_INSTALL_RPATH "${TTK_BINARY_INSTALL_DIR}/../")

  # add the common package
  ttk_add_baseCode_package(common)

  if("${TTK_BUILD_MODE}" MATCHES "VTK-GUI")
    set(withVTK true)
    ttk_add_vtkWrapper_package(ttkWrapper)
  endif("${TTK_BUILD_MODE}" MATCHES "VTK-GUI")

  if("${TTK_BUILD_MODE}" MATCHES "VTK-CMD")
    set(withVTK true)
    ttk_add_vtkWrapper_package(ttkWrapper)
  endif("${TTK_BUILD_MODE}" MATCHES "VTK-CMD")

  if(NOT "${TTK_BUILD_MODE}" MATCHES "ParaView")
    if(withVTK)
      if(NOT VTK_USE_FILE)
        find_package(VTK REQUIRED)
        include(${VTK_USE_FILE})  
        set(PROJECT_DEP "${PROJECT_DEP} ${VTK_LIBRARIES}" 
          CACHE INTERNAL "PROJECT_DEP")
        set(PROJECT_FLAGS "${PROJECT_FLAGS} -DwithVTK"
          CACHE INTERNAL "PROJECT_FLAGS")

        if("${TTK_BUILD_MODE}" MATCHES "VTK-CMD")
          ttk_add_vtkWrapper_package(ttkProgramBase)
        endif("${TTK_BUILD_MODE}" MATCHES "VTK-CMD")
          
        if("${TTK_BUILD_MODE}" MATCHES "VTK-GUI")
          ttk_add_vtkWrapper_package(ttkUserInterfaceBase)
          ttk_add_vtkWrapper_package(ttkTextureMapFromField)
          ttk_add_vtkWrapper_package(ttkWRLExporter)
        endif("${TTK_BUILD_MODE}" MATCHES "VTK-GUI")
      endif(NOT VTK_USE_FILE)
    endif(withVTK)
  endif(NOT "${TTK_BUILD_MODE}" MATCHES "ParaView")

  if("${TTK_BUILD_MODE}" MATCHES "ParaView")
    find_package(ParaView REQUIRED)
    include(${PARAVIEW_USE_FILE})
    set(TTK_PLUGIN_INSTALL_DIR "${VTK_INSTALL_PREFIX}/lib/paraview-${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}/plugins"
      CACHE PATH
      "Directory where TTK ParaView plugins will be installed")
    set(withVTK true)
    ttk_add_vtkWrapper_package(ttkWrapper)
  endif("${TTK_BUILD_MODE}" MATCHES "ParaView")

  # messages on dependencies
  message(STATUS
"ttk -------------------------------------------------------------------------")
  message(STATUS "baseCode path: ${BASECODE_DIR}")

  # optional dependency on VTK
  message(STATUS "Binary installation path: ${TTK_BINARY_INSTALL_DIR}")
  message(STATUS "Plugin installation path: ${TTK_PLUGIN_INSTALL_DIR}")
  message(STATUS
"-----------------------------------------------------------------------------")

  set(TTK_PKG true)

endif(NOT TTK_PKG)
