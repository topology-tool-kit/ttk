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
#     LINK <libraries to link>: DEPRECATED)
#
function(ttk_add_paraview_plugin library)
	cmake_parse_arguments(ARG "" "PLUGIN_XML;PLUGIN_VERSION"
		"SOURCES;LINK" ${ARGN})

	if(DEFINED ARG_PLUGIN_VERSION)
		message(DEPRECATED ": VERSION is deprecated in ttk_add_paraview_plugin (${CMAKE_CURRENT_SOURCE_DIR}) ")
	endif()

	set(TTK_PV_XML     "${CMAKE_CURRENT_SOURCE_DIR}/${ARG_PLUGIN_XML};${TTK_PV_XML}" CACHE INTERNAL "")
	set(TTK_PV_SOURCES "${ARG_SOURCES};${TTK_PV_SOURCES}"                            CACHE INTERNAL "")
	set(TTK_PV_LINKS   "${ARG_LINK};${TTK_PV_LINKS}"                                CACHE INTERNAL "")
	set(TTK_PV_DIRS    "${CMAKE_CURRENT_SOURCE_DIR};${TTK_PV_DIRS}"                  CACHE INTERNAL "")
endfunction()
