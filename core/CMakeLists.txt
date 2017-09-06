find_package(Doxygen)
if (DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ttk.doxygen
		${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen)
	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif()

add_subdirectory(baseCode)
add_subdirectory(vtkWrappers)

