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
	add_library(${library} ${ARG_SOURCES})
	target_link_libraries(${library} PUBLIC ${VTK_LIBRARIES} ${ARG_LINK})

	target_include_directories(${library} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:include/ttk/vtk>)

	install(TARGETS ${library}
		EXPORT TTKVTKConfig
		RUNTIME DESTINATION bin/ttk/vtk/
		ARCHIVE DESTINATION lib/ttk/vtk
		LIBRARY DESTINATION lib/ttk/vtk)

	install(FILES ${ARG_HEADERS} DESTINATION include/ttk/vtk)

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

		install(TARGETS ${plugin_name} DESTINATION plugin)
	endif()
endfunction()

add_subdirectory(ttkBlank)
add_subdirectory(ttkBottleneckDistance)
add_subdirectory(ttkCellDataConverter)
add_subdirectory(ttkComponentSize)
add_subdirectory(ttkContinuousScatterPlot)
add_subdirectory(ttkContourForests)
add_subdirectory(ttkDistanceField)
add_subdirectory(ttkFiber)
add_subdirectory(ttkFiberSurface)
add_subdirectory(ttkGeometrySmoother)
add_subdirectory(ttkIdentifierRandomizer)
add_subdirectory(ttkIdentifiers)
add_subdirectory(ttkIntegralLines)
add_subdirectory(ttkJacobiSet)
add_subdirectory(ttkLDistance)
add_subdirectory(ttkMandatoryCriticalPoints)
add_subdirectory(ttkMeshSubdivision)
add_subdirectory(ttkMorseSmaleComplex)
add_subdirectory(ttkPersistenceCurve)
add_subdirectory(ttkPersistenceDiagram)
add_subdirectory(ttkPointDataConverter)
add_subdirectory(ttkProgramBase)
add_subdirectory(ttkProjectionFromField)
add_subdirectory(ttkRangePolygon)
add_subdirectory(ttkReebSpace)
add_subdirectory(ttkScalarFieldCriticalPoints)
add_subdirectory(ttkScalarFieldNormalizer)
add_subdirectory(ttkScalarFieldSmoother)
add_subdirectory(ttkSphereFromPoint)
add_subdirectory(ttkTextureMapFromField)
add_subdirectory(ttkTopologicalSimplification)
add_subdirectory(ttkTriangulation)
add_subdirectory(ttkTriangulationRequest)
add_subdirectory(ttkUncertainDataEstimator)
add_subdirectory(ttkUserInterfaceBase)
add_subdirectory(ttkWRLExporter)

install(EXPORT TTKVTKConfig
	NAMESPACE ttk::vtk::
	DESTINATION lib/cmake/ttk)

