add_subdirectory(base)

if(TTK_BUILD_DOCUMENTATION)
  if(DOXYGEN_FOUND)
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
    install(
      DIRECTORY
        ${CMAKE_SOURCE_DIR}/doc/img
      DESTINATION
        ${CMAKE_INSTALL_DATAROOTDIR}/doc/ttk
        )
  endif()
endif()

