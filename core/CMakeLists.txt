# Add compile flags and defintions to the target
# according to the options selected by the user.
#
# Usage:
# ttk_set_compile_options(<library_name>)
#
function(ttk_set_compile_options library)
	target_compile_options(${library} PRIVATE -Wall)

	if (TTK_ENABLE_KAMIKAZE)
		target_compile_definitions(${library} PUBLIC TTK_ENABLE_KAMIKAZE)
	endif()

	if (TTK_ENABLE_CPU_OPTIMIZATION)
		target_compile_options(${library} PRIVATE -march=native)
	endif()

	if (TTK_ENABLE_OPENMP)
		target_compile_definitions(${library} PUBLIC TTK_ENABLE_OPENMP)
		target_compile_options(${library} PUBLIC ${OpenMP_CXX_FLAGS})
		target_link_libraries(${library} PUBLIC ${OpenMP_CXX_FLAGS})
	endif()

	if (TTK_ENABLE_MPI)
		target_compile_definitions(${library} PUBLIC TTK_ENABLE_MPI)
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

if (NOT APPLE)
	option(TTK_ENABLE_OPENMP "Enable OpenMP support" TRUE)
endif()
option(TTK_ENABLE_MPI "Enable MPI support" FALSE)

if (TTK_ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
endif()

if (TTK_ENABLE_MPI)
	find_package(MPI REQUIRED)
endif()

add_subdirectory(base)
if (TTK_BUILD_VTK_WRAPPERS)
	add_subdirectory(vtk)
endif()

