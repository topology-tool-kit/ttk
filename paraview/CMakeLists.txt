# Fill TTK_MODULES, TTK_MODULES_FILES
# These lists are then used by singlePlugin/CMakeLists.txt
# The goal is to create a single plugin for parview
set(TTK_XMLS             "" CACHE INTERNAL "")
set(TTK_MODULES          "" CACHE INTERNAL "")
file(GLOB XML_FILES xmls/*.xml)

foreach(XML_FILE ${XML_FILES})
  get_filename_component(MODULE_NAME ${XML_FILE} NAME_WLE)
  ttk_register_pv_filter(ttk${MODULE_NAME} xmls/${MODULE_NAME}.xml)
endforeach()

set("PARAVIEW_PLUGIN_ENABLE_TopologyToolKit"
  ON
  CACHE
  BOOL
  "Build TTK TopologyToolKit filter"
  FORCE
  )

paraview_plugin_scan(
  PLUGIN_FILES
    "singlePlugin/paraview.plugin"
  PROVIDES_PLUGINS
    TTK_PLUGINS
  REQUIRES_MODULES
    TTK_REQUIRED
  )

foreach(REQ IN LISTS TTK_REQUIRED)
  if(NOT TARGET "${REQ}")
    message(FATAL_ERROR "Missing required module: ${REQ}")
  endif()
endforeach()

paraview_plugin_build(
  TARGET
    TTKPVTargets
  PLUGINS
    ${TTK_PLUGINS}
  AUTOLOAD
    ${TTK_PLUGINS}
  INSTALL_HEADERS
    OFF
  RUNTIME_DESTINATION
    "${CMAKE_INSTALL_BINDIR}/plugins"
  LIBRARY_DESTINATION
    "${CMAKE_INSTALL_BINDIR}/plugins"
  )

install(
  TARGETS
    TTKPVTargets
  DESTINATION
    "${TTK_INSTALL_PLUGIN_DIR}"
  )