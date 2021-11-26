set(CPACK_PACKAGE_NAME "TTK")
set(CPACK_PACKAGE_FILE_NAME "ttk")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Topology ToolKit")
set(CPACK_PACKAGE_VERSION_MAJOR ${TTK_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${TTK_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${TTK_VERSION_PATCH})
set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_CONTACT "Julien Tierny <julien.tierny@sorbonne-universite.fr>")
set(CPACK_PACKAGE_VENDOR "CNRS, Sorbonne University and contributors")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://topology-tool-kit.github.io/")
if(NOT APPLE)
  set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)
  set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
else()
  # macOS needs license & readme files ending with .txt
  configure_file(LICENSE License.txt COPYONLY)
  configure_file(README.md Readme.txt COPYONLY)
  set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_BINARY_DIR}/License.txt)
  set(CPACK_RESOURCE_FILE_README ${PROJECT_BINARY_DIR}/Readme.txt)
endif()
set(CPACK_DEBIAN_PACKAGE_DEPENDS
  "ttk-paraview (= 5.9.1), python3-sklearn, libboost-system-dev, python3-dev, libgraphviz-dev, libsqlite3-dev, libgl1-mesa-dev")
# autogenerate dependency information
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
# package will be installed under %ProgramFiles%\${CPACK_PACKAGE_INSTALL_DIRECTORY} on Windows
set(CPACK_PACKAGE_INSTALL_DIRECTORY "TTK")
# let the installer uninstall previous installations on Windows
set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
# generate components, fix productbuild packaging for macOS
if(APPLE)
  set(CPACK_COMPONENTS_ALL Unspecified python development)
endif()
# embed Visual C++ redistribuable for Windows + set environment variables
if(WIN32 AND EXISTS "${CMAKE_SOURCE_DIR}/vc_redist.x64.exe")
  install(FILES "${CMAKE_SOURCE_DIR}/scripts/set_windows_env_vars.ps1" DESTINATION bin)
  install(PROGRAMS "${CMAKE_SOURCE_DIR}/vc_redist.x64.exe" DESTINATION bin)
  set(CPACK_VERBATIM_VARIABLES TRUE)
  list(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS
    " ExecWait 'powershell -ExecutionPolicy Bypass -File \"$INSTDIR\\bin\\set_windows_env_vars.ps1\" install'
      ExecWait '$INSTDIR\\bin\\vc_redist.x64.exe /passive /norestart'
      Delete '$INSTDIR\\bin\\vc_redist.x64.exe' ")
  list(APPEND CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS
    " ExecWait 'powershell -ExecutionPolicy Bypass -File \"$INSTDIR\\bin\\set_windows_env_vars.ps1\"' ")
endif()
include(CPack)
