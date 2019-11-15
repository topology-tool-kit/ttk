# fix for typo in ParaView 5.7
# https://discourse.paraview.org/t/install-headers-for-vtk-modules-build-by-paraview-5-7-0/2721
set(_paraview_add_plugin_MODULES_FILES TRUE)

paraview_add_plugin(TopologyToolKit
  VERSION "1.0"
  REQUIRED_ON_CLIENT
  REQUIRED_ON_SERVER
  MODULES
    ${TTK_MODULES}
  MODULE_FILES
    ${TTK_VTK_MODULE_FILES}
  MODULE_ARGS # BUG: see above
    INSTALL_HEADERS ON
  SERVER_MANAGER_XML
    ${TTK_XMLS}
  )

if(NOT "${TTK_INSTALL_PLUGIN_DIR}" STREQUAL "")
  install(
    TARGETS
      TopologyToolKit
    DESTINATION
      "${TTK_INSTALL_PLUGIN_DIR}"
    )
endif()
