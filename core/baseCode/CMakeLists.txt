# Function to create the library and install targets
# for a TTK BaseCode library.
#
# Usage:
# ttk_add_base_library(<library_name> <source list> HEADERS <headers to install>)
#
function(ttk_add_base_library library)
	set(SOURCES "")
	set(HEADER_LIST "")

	set(CURRENT_LIST SOURCES)
	foreach (arg ${ARGN})
		if ("${arg}" STREQUAL "HEADERS")
			set(CURRENT_LIST HEADER_LIST)
		else()
			list(APPEND ${CURRENT_LIST} ${arg})
		endif()
	endforeach()

	add_library(${library} SHARED ${SOURCES})
	target_include_directories(${library} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:include/ttk/base>)

	install(TARGETS ${library}
		EXPORT TTKBaseConfig
		RUNTIME DESTINATION bin/ttk/base/
		ARCHIVE DESTINATION lib/ttk/base
		LIBRARY DESTINATION lib/ttk/base)

	# Generate the ttk::base namespaced target name
	string(SUBSTRING "${library}" 0 1 first_letter)
	string(SUBSTRING "${library}" 1 -1 rest_of_name)
	string(TOUPPER "${first_letter}" first_letter)
	string(CONCAT namespaced_library ${first_letter} ${rest_of_name})
	add_library(ttk::base::${namespaced_library} ALIAS ${library})

	install(FILES ${HEADER_LIST} DESTINATION include/ttk/base)
endfunction()

# Function to create a header only template library and install targets
# for a TTK BaseCode template library.
#
# Usage:
# ttk_add_base_template_library(<library_name> <headers list>)
#
function(ttk_add_base_template_library library)
	set(SOURCES ${ARGN})

	add_library(${library} INTERFACE)
	target_include_directories(${library} INTERFACE
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:include/ttk/base>)

	install(TARGETS ${library}
		EXPORT TTKBaseConfig
		RUNTIME DESTINATION bin/ttk/base/
		ARCHIVE DESTINATION lib/ttk/base
		LIBRARY DESTINATION lib/ttk/base)

	# Generate the ttk::base namespaced target name
	string(SUBSTRING "${library}" 0 1 first_letter)
	string(SUBSTRING "${library}" 1 -1 rest_of_name)
	string(TOUPPER "${first_letter}" first_letter)
	string(CONCAT namespaced_library ${first_letter} ${rest_of_name})
	add_library(ttk::base::${namespaced_library} ALIAS ${library})

	install(FILES ${ARGN} DESTINATION include/ttk/base)
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

