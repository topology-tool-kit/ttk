include(functions.cmake)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
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
target_link_libraries(ttkAll INTERFACE ${VTKWRAPPER_LIB_LIST})

install(TARGETS ttkAll
	EXPORT TTKVTKTargets
	RUNTIME DESTINATION bin/ttk
	ARCHIVE DESTINATION lib/ttk
	LIBRARY DESTINATION lib/ttk)

install(EXPORT TTKVTKTargets
	NAMESPACE ttk::vtk::
	DESTINATION lib/cmake/ttk)

configure_file(TTKVTKConfig.cmake.in TTKVTKConfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/TTKVTKConfig.cmake"
	DESTINATION lib/cmake/ttk)

