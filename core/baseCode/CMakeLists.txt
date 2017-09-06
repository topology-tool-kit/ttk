# Function to create the library and install targets
# for a TTK BaseCode library.
#
# Usage:
# ttk_add_base_library(<library_name>
#     SOURCES <source list>
#     HEADERS <headers to install>
#     LINK <libraries to link>)
#
function(ttk_add_base_library library)
	cmake_parse_arguments(ARG "" "" "SOURCES;HEADERS;LINK" ${ARGN})

	add_library(${library} ${ARG_SOURCES})
	target_include_directories(${library} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:include/ttk/base>)
	if (ARG_LINK)
		target_link_libraries(${library} PUBLIC ${ARG_LINK})
	endif()

	install(TARGETS ${library}
		EXPORT TTKBaseConfig
		RUNTIME DESTINATION bin/ttk/base/
		ARCHIVE DESTINATION lib/ttk/base
		LIBRARY DESTINATION lib/ttk/base)

	if (ARG_HEADERS)
		install(FILES ${ARG_HEADERS} DESTINATION include/ttk/base)
	endif()
endfunction()

# Function to create a header only template library and install targets
# for a TTK BaseCode template library. Note that here "libraries to link"
# has a bit of a different meaning. Since it's a template library no code
# is compiled or linked, but this "linking" is used to express dependencies
# on other libraries that users of this library will also need.
#
# Usage:
# ttk_add_base_template_library(<library_name>
#     SOURCES <headers list>
#     LINK <libraries to link>)
#
function(ttk_add_base_template_library library)
	cmake_parse_arguments(ARG "" "" "SOURCES;LINK" ${ARGN})

	add_library(${library} INTERFACE)
	target_include_directories(${library} INTERFACE
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:include/ttk/base>)

	install(TARGETS ${library}
		EXPORT TTKBaseConfig
		RUNTIME DESTINATION bin/ttk/base/
		ARCHIVE DESTINATION lib/ttk/base
		LIBRARY DESTINATION lib/ttk/base)
	if (ARG_LINK)
		target_link_libraries(${library} INTERFACE ${ARG_LINK})
	endif()

	install(FILES ${ARG_SOURCES} DESTINATION include/ttk/base)
endfunction()

add_subdirectory(abstractMorseSmaleComplex)
add_subdirectory(abstractTriangulation)
add_subdirectory(blank)
add_subdirectory(bottleneckDistance)
add_subdirectory(common)
add_subdirectory(continuousScatterPlot)
add_subdirectory(contourForests)
add_subdirectory(contourForestsTree)
add_subdirectory(contourTree)
add_subdirectory(discreteGradient)
add_subdirectory(distanceField)
add_subdirectory(explicitTriangulation)
add_subdirectory(fiberSurface)
add_subdirectory(geometry)
add_subdirectory(implicitTriangulation)
add_subdirectory(integralLines)
add_subdirectory(jacobiSet)
add_subdirectory(lDistance)
add_subdirectory(lowestCommonAncestor)
add_subdirectory(mandatoryCriticalPoints)
add_subdirectory(morseSmaleComplex)
add_subdirectory(morseSmaleComplex2D)
add_subdirectory(morseSmaleComplex3D)
add_subdirectory(persistenceCurve)
add_subdirectory(persistenceDiagram)
add_subdirectory(rangeDrivenOctree)
add_subdirectory(rangeMinimumQuery)
add_subdirectory(reebSpace)
add_subdirectory(scalarFieldCriticalPoints)
add_subdirectory(scalarFieldSmoother)
add_subdirectory(skeleton)
add_subdirectory(topologicalSimplification)
add_subdirectory(triangulation)
add_subdirectory(uncertainDataEstimator)
add_subdirectory(unionFind)

install(EXPORT TTKBaseConfig
	NAMESPACE ttk::base::
	DESTINATION lib/cmake/ttk)

