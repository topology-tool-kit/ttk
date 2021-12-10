add_subdirectory(base)

if(TTK_BUILD_DOCUMENTATION)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_INCLUDE_DIRS)
    if(TARGET VTK::CommonCore)
      get_target_property(DOXYGEN_INCLUDE_DIRS VTK::CommonCore INTERFACE_INCLUDE_DIRECTORIES)
    endif()
    # to find welcomeMsg.inl and welcomeLogo.inl
    list(APPEND DOXYGEN_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/core/base/common")
    # Doxygen needs space-separated paths
    string (REPLACE ";" " " DOXYGEN_INCLUDE_DIRS "${DOXYGEN_INCLUDE_DIRS}")
    # to find splash.png
    set(DOXYGEN_IMAGE_PATH "${PROJECT_SOURCE_DIR}/doc/img")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ttk.doxygen
      ${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen)
    add_custom_target(doc
      ALL
        ${DOXYGEN_EXECUTABLE}
        ${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen
      WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT
        "Generating API documentation with Doxygen"
      VERBATIM
      )
    install(
      DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/doc/html
      DESTINATION
        ${CMAKE_INSTALL_DATAROOTDIR}/doc/ttk
        )
  endif()
endif()
