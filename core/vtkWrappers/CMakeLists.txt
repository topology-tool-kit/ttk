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
	set(HEADERS "")

	set(CURRENT_LIST SOURCES)
	foreach (arg ${ARGN})
		if ("${arg}" STREQUAL "HEADERS")
			set(CURRENT_LIST HEADERS)
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

	install(FILES ${HEADERS} DESTINATION include/ttk/vtk)
endfunction()

add_subdirectory(ttkWrapper)
add_subdirectory(vtkBlank)
add_subdirectory(vtkCellDataConverter)
add_subdirectory(vtkComponentSize)
add_subdirectory(vtkContinuousScatterPlot)
add_subdirectory(vtkContourForests)
add_subdirectory(vtkDistanceField)
add_subdirectory(vtkFiber)
add_subdirectory(vtkFiberSurface)
add_subdirectory(vtkGeometrySmoother)
add_subdirectory(vtkIdentifierRandomizer)
add_subdirectory(vtkIdentifiers)
add_subdirectory(vtkIntegralLines)
add_subdirectory(vtkJacobiSet)
add_subdirectory(vtkMandatoryCriticalPoints)
add_subdirectory(vtkMeshSubdivision)
add_subdirectory(vtkMorseSmaleComplex)
add_subdirectory(vtkPersistenceCurve)
add_subdirectory(vtkPersistenceDiagram)
add_subdirectory(vtkPointDataConverter)
# TODO: Is program base used? It was commented out
add_subdirectory(vtkProgramBase)
add_subdirectory(vtkProjectionFromField)
add_subdirectory(vtkRangePolygon)
add_subdirectory(vtkReebSpace)
add_subdirectory(vtkScalarFieldCriticalPoints)
add_subdirectory(vtkScalarFieldNormalizer)
add_subdirectory(vtkScalarFieldSmoother)
add_subdirectory(vtkSphereFromPoint)
add_subdirectory(vtkTextureMapFromField)
add_subdirectory(vtkTopologicalSimplification)
add_subdirectory(vtkTriangulation)
add_subdirectory(vtkTriangulationRequest)
add_subdirectory(vtkUncertainDataEstimator)
add_subdirectory(vtkUserInterfaceBase)
add_subdirectory(vtkWRLExporter)


install(EXPORT TTKVTKConfig
	NAMESPACE ttk::vtk::
	DESTINATION lib/cmake/ttk)

