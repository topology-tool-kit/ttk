find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

# Function to create the library and install targets
# for a TTK VTK Wrapper library.
#
# Usage:
# ttk_add_vtk_library(<library_name> <source list> HEADERS <headers to install>)
#
function(ttk_add_vtk_library library)
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

	add_library(${library} ${SOURCES})
	target_link_libraries(${library} ${VTK_LIBRARIES})
	target_include_directories(${library} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:include/ttk/vtk>)

	install(TARGETS ${library}
		EXPORT TTKVTKConfig
		RUNTIME DESTINATION bin/ttk/vtk/
		ARCHIVE DESTINATION lib/ttk/vtk
		LIBRARY DESTINATION lib/ttk/vtk)

	# Generate the ttk::vtk namespaced target name
	string(SUBSTRING "${library}" 3 -1 unprefixed_library)
	add_library(ttk::vtk::${unprefixed_library} ALIAS ${library})

	install(FILES ${HEADER_LIST} DESTINATION include/ttk/vtk)
endfunction()

add_subdirectory(ttkWrapper)
add_subdirectory(ttkBlank)
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
add_subdirectory(ttkMandatoryCriticalPoints)
add_subdirectory(ttkMeshSubdivision)
add_subdirectory(ttkMorseSmaleComplex)
add_subdirectory(ttkPersistenceCurve)
add_subdirectory(ttkPersistenceDiagram)
add_subdirectory(ttkPointDataConverter)
# TODO: Is program base used? It was commented out
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

