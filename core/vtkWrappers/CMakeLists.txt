# Function to create the library and install targets
# for a TTK VTK Wrapper library. The VTK Libraries will
# automatically be linked.
# Options:
# SOURCES: Specify the list of source files for the library
# HEADERS: Specify the list of header files to install for the library
# LINK: Specify the link dependencies of the library
# PLUGIN_VERSION: Specify ParaView plugin version
# PLUGIN_XML: Specify ParaView plugin XML file
#
# To build a VTK Wrapper library:
# ttk_add_vtk_library(<library_name>
#     SOURCES <source list>
#     HEADERS <headers to install>
#     LINK <libraries to link>)
#
# To build a VTK Wrapper library and optional ParaView plugin:
# (the plugin will only be built if TTK_BUILD_PARAVIEW_PLUGINS is enabled)
#
# ttk_add_vtk_library(<library_name>
#     SOURCES <source list>
#     HEADERS <headers to install>
#     PLUGIN_VERSION <version>
#     PLUGIN_XML <XML file for the plugin>
#     LINK <libraries to link>)
#
function(ttk_add_vtk_library library)
	cmake_parse_arguments(ARG "" "PLUGIN_XML;PLUGIN_VERSION"
		"SOURCES;HEADERS;LINK" ${ARGN})

	# Build the VTK Wrapper
	add_library(${library} SHARED ${ARG_SOURCES})
	target_link_libraries(${library} PUBLIC ${VTK_LIBRARIES} ${ARG_LINK})

	target_include_directories(${library} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:include/ttk/vtk>)

	ttk_set_compile_options(${library})

	install(TARGETS ${library}
		EXPORT TTKVTKConfig
		RUNTIME DESTINATION bin/ttk
		ARCHIVE DESTINATION lib/ttk
		LIBRARY DESTINATION lib/ttk)

	if (ARG_HEADERS)
		install(FILES ${ARG_HEADERS} DESTINATION include/ttk/vtk)
	endif()

	# Build the ParaView plugin, if we're building with ParaView
	# and this wrapper library has a plugin
	if (TTK_BUILD_PARAVIEW_PLUGINS AND ARG_PLUGIN_XML)
		set(plugin_name "${library}Plugin")
		add_paraview_plugin(${plugin_name} ${ARG_PLUGIN_VERSION}
			SERVER_MANAGER_XML ${ARG_PLUGIN_XML}
			SERVER_MANAGER_SOURCES ${ARG_SOURCES})
		target_link_libraries(${plugin_name} PUBLIC ${VTK_LIBRARIES} ${ARG_LINK})

		target_include_directories(${plugin_name} PUBLIC
			$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

		install(TARGETS ${plugin_name} DESTINATION ${TTK_PLUGIN_INSTALL_DIR})
	endif()
endfunction()

option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" OFF)
if (NOT TTK_BUILD_PARAVIEW_PLUGINS)
	option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" OFF)
	find_package(VTK REQUIRED)
	include(${VTK_USE_FILE})
	include_directories(${VTK_INCLUDE_DIRS})
else()
	set(TTK_BUILD_VTK_WRAPPERS ON CACHE BOOL "Build the TTK VTK Wrappers" FORCE)
	find_package(ParaView REQUIRED)
	include(${PARAVIEW_USE_FILE})
	include_directories(${PARAVIEW_INCLUDE_DIRS})
	set(TTK_PLUGIN_INSTALL_DIR
		"${VTK_INSTALL_PREFIX}/lib/paraview-${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}/plugins"
		CACHE PATH
		"Directory where TTK ParaView plugins will be installed")
endif()

file(GLOB VTKWRAPPER_DIRS ttk*)
foreach(VTKWRAPPER ${VTKWRAPPER_DIRS})
	if (IS_DIRECTORY ${VTKWRAPPER})
		add_subdirectory(${VTKWRAPPER})
	endif()
endforeach()

install(EXPORT TTKVTKConfig
	NAMESPACE ttk::vtk::
	DESTINATION lib/cmake/ttk)

