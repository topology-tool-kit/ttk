# Function to create the ParaView plugin for a TTK ParaView plugin

# Options:
# SOURCES: Specify the list of source files for the library
# LINK: Specify the link dependencies of the library
# PLUGIN_VERSION: Specify ParaView plugin version
# PLUGIN_XML: Specify ParaView plugin XML file
#
# ttk_add_paraview_plugin(<library_name>
#     SOURCES <source list>
#     PLUGIN_VERSION <version>
#     PLUGIN_XML <XML file for the plugin>
#     LINK <libraries to link>)
#
function(ttk_add_paraview_plugin library)
	cmake_parse_arguments(ARG "" "PLUGIN_XML;PLUGIN_VERSION"
		"SOURCES;LINK" ${ARGN})

	# Build the ParaView plugin, if we're building with ParaView
	# and this wrapper library has a plugin
	set(plugin_name "${library}Plugin")
	add_paraview_plugin(${plugin_name} ${ARG_PLUGIN_VERSION}
		SERVER_MANAGER_XML ${ARG_PLUGIN_XML}
		SERVER_MANAGER_SOURCES ${ARG_SOURCES} ${TTK_TRIANGULATION_SRCS})
	target_link_libraries(${plugin_name} PUBLIC ${VTK_LIBRARIES} ${ARG_LINK} ttkPVTriangulation)

	if (MSVC)
		target_compile_definitions(${plugin_name} PUBLIC TTK_PLUGIN)
	endif()

	target_include_directories(${plugin_name} PUBLIC
		$<BUILD_INTERFACE:${VTKWRAPPER_DIR}/${library}>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

	install(TARGETS ${plugin_name} DESTINATION ${TTK_INSTALL_PLUGIN_DIR})
endfunction()
