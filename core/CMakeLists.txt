# Add compile flags and defintions to the target
# according to the options selected by the user.
#
# Usage:
# ttk_set_compile_options(<library_name>)
#
function(ttk_set_compile_options library)
	target_compile_options(${library} PRIVATE -Wall)

	if (TTK_WITH_KAMIKAZE)
		target_compile_definitions(${library} PUBLIC TTK_WITH_KAMIKAZE)
	endif()

	if (TTK_ENABLE_CPU_OPTIMIZATION)
		target_compile_options(${library} PRIVATE -march=native)
	endif()

	if (TTK_WITH_OPENMP)
		target_compile_definitions(${library} PUBLIC TTK_WITH_OPENMP)
		target_compile_options(${library} PUBLIC ${OpenMP_CXX_FLAGS})
		target_link_libraries(${library} PUBLIC ${OpenMP_CXX_FLAGS})
	endif()

	if (TTK_WITH_MPI)
		target_compile_definitions(${library} PUBLIC TTK_WITH_MPI)
		target_include_directories(${library} PUBLIC ${MPI_CXX_INCLUDE_PATH})
		target_link_libraries(${library} PUBLIC ${MPI_CXX_LIBRARIES})
	endif()
endfunction()

find_package(Doxygen)
if (DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ttk.doxygen
		${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen)
	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif()

option(TTK_WITH_KAMIKAZE "Enable Kamikaze compilation mode" FALSE)
option(TTK_ENABLE_CPU_OPTIMIZATION "Enable native CPU optimizations" TRUE)
option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" ON)
option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" ON)
if (NOT APPLE)
	option(TTK_WITH_OPENMP "Enable OpenMP support" TRUE)
endif()
option(TTK_WITH_MPI "Enable MPI support" FALSE)

if (TTK_WITH_OPENMP)
	find_package(OpenMP REQUIRED)
endif()

if (TTK_WITH_MPI)
	find_package(MPI REQUIRED)
endif()

if (TTK_BUILD_PARAVIEW_PLUGINS)
	set(TTK_BUILD_VTK_WRAPPERS ON CACHE BOOL "Build the TTK VTK Wrappers" FORCE)
endif()

add_subdirectory(baseCode)
if (TTK_BUILD_VTK_WRAPPERS)
	add_subdirectory(vtkWrappers)
endif()

