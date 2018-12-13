include(functions.cmake)

find_package(ParaView REQUIRED)
include(${PARAVIEW_USE_FILE})
include_directories(${PARAVIEW_INCLUDE_DIRS})

if (MSVC)
set(TTK_INSTALL_PLUGIN_DIR
	"$ENV{PARAVIEW_BIN}/plugins"
	CACHE PATH
	"Directory where TTK ParaView plugins will be installed")
else()
set(TTK_INSTALL_PLUGIN_DIR
	"${VTK_INSTALL_PREFIX}/lib/plugins"
	CACHE PATH
	"Directory where TTK ParaView plugins will be installed")
endif()

set(VTKWRAPPER_DIR "${CMAKE_CURRENT_LIST_DIR}/../core/vtk/")

# Build the ttkTriangulation library against the ParaView VTK
add_library(ttkPVTriangulation STATIC
	${VTKWRAPPER_DIR}/ttkTriangulation/ttkTriangulation.cpp
	${VTKWRAPPER_DIR}/ttkTriangulation/ttkWrapper.cpp)
target_link_libraries(ttkPVTriangulation PUBLIC ${VTK_LIBRARIES} triangulation)
set_target_properties(ttkPVTriangulation PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_include_directories(ttkPVTriangulation PUBLIC
	$<BUILD_INTERFACE:${VTKWRAPPER_DIR}/ttkTriangulation>)
if (MSVC)
	target_compile_definitions(ttkPVTriangulation PUBLIC TTK_PLUGIN)
endif()

# traverse sub-directories
file(GLOB PARAVIEW_PLUGIN_DIRS *)
foreach(PARAVIEW_PLUGIN ${PARAVIEW_PLUGIN_DIRS})
	if (IS_DIRECTORY ${PARAVIEW_PLUGIN})
		if (NOT ${PARAVIEW_PLUGIN} MATCHES "patch" )
			add_subdirectory(${PARAVIEW_PLUGIN})
		endif()
	endif()
endforeach()


