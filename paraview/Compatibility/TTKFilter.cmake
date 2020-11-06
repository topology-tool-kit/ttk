option(PARAVIEW_PLUGIN_ENABLE_Compatibility "Add filter definitions for backwards compatibility with v0.9.8.9 to the TopologyToolkit plugin." OFF)

if(PARAVIEW_PLUGIN_ENABLE_Compatibility)
  list(APPEND TTK_XMLS "${CMAKE_CURRENT_LIST_DIR}/Compatibility.xml")
endif()