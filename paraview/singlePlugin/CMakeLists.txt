# fix for typo in ParaView 5.7
# https://discourse.paraview.org/t/install-headers-for-vtk-modules-build-by-paraview-5-7-0/2721
set(_paraview_add_plugin_MODULES_FILES TRUE)

list(REMOVE_DUPLICATES TTK_MODULES)

paraview_add_plugin(TopologyToolKit
  VERSION ${CMAKE_PROJECT_VERSION}
  REQUIRED_ON_CLIENT
  REQUIRED_ON_SERVER
  MODULES
    ${TTK_MODULES}
  SERVER_MANAGER_XML
    ${TTK_XMLS}
  )

if(NOT "${TTK_INSTALL_PLUGIN_DIR}" STREQUAL "")
  ttk_set_paraview_rpath(TopologyToolKit)

  install(
    TARGETS
      TopologyToolKit
    DESTINATION
      "${TTK_INSTALL_PLUGIN_DIR}/${TTK_PLUGIN_SUBDIR}"
    )
endif()
