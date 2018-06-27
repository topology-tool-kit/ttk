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

  if(NOT MSVC)
    if (TTK_ENABLE_CPU_OPTIMIZATION)
      target_compile_options(${library} 
        PRIVATE $<$<CONFIG:Release>:-march=native -O3 -Wfatal-errors>)
    endif()

    target_compile_options(${library} PRIVATE $<$<CONFIG:Debug>:-O0 -g -pg>)
  endif()

	if (TTK_ENABLE_OPENMP)
		target_compile_definitions(${library} PUBLIC TTK_ENABLE_OPENMP)
		target_compile_options(${library} PUBLIC ${OpenMP_CXX_FLAGS})
		target_link_libraries(${library} PUBLIC ${OpenMP_CXX_LIBRARIES})
	endif()

	if (TTK_ENABLE_MPI)
		target_compile_definitions(${library} PUBLIC TTK_ENABLE_MPI)
		target_include_directories(${library} PUBLIC ${MPI_CXX_INCLUDE_PATH})
		target_link_libraries(${library} PUBLIC ${MPI_CXX_LIBRARIES})
	endif()

	if (TTK_ENABLE_ZFP)
		target_compile_definitions(${library} PUBLIC TTK_ENABLE_ZFP)
		target_include_directories(${library} PUBLIC ${ZFP_INCLUDE_DIR})
		target_link_libraries(${library} PUBLIC ${ZFP_LIBRARY})
	endif()
	
	if (TTK_ENABLE_ZLIB)
		target_compile_definitions(${library} PUBLIC TTK_ENABLE_ZLIB)
		target_include_directories(${library} PUBLIC ${ZLIB_INCLUDE_DIR})
		target_link_libraries(${library} PUBLIC ${ZLIB_LIBRARY})
	endif()

endfunction()

function(ttk_print_summary)
	message(STATUS "ttk options -----------------------------------------------------------------")
	message(STATUS "TTK_ENABLE_CPU_OPTIMIZATION: ${TTK_ENABLE_CPU_OPTIMIZATION}")
	message(STATUS "TTK_ENABLE_KAMIKAZE: ${TTK_ENABLE_KAMIKAZE}")
	message(STATUS "TTK_ENABLE_OPENMP: ${TTK_ENABLE_OPENMP}")
	message(STATUS "TTK_ENABLE_MPI: ${TTK_ENABLE_MPI}")
  message(STATUS "TTK_ENABLE_ZFP: ${TTK_ENABLE_ZFP}")
  message(STATUS "TTK_ENABLE_ZLIB: ${TTK_ENABLE_ZLIB}")
	message(STATUS "ttk build -------------------------------------------------------------------")
	message(STATUS "TTK_BUILD_DOCUMENTATION: ${TTK_BUILD_DOCUMENTATION}")
	if(TTK_BUILD_DOCUMENTATION)
		message(STATUS "  DOXYGEN_EXECUTABLE: ${DOXYGEN_EXECUTABLE}")
	endif()
	message(STATUS "TTK_BUILD_PARAVIEW_PLUGINS: ${TTK_BUILD_PARAVIEW_PLUGINS}")
	if(TTK_BUILD_PARAVIEW_PLUGINS)
		message(STATUS "  ParaView_DIR: ${ParaView_DIR}")
	endif()
	message(STATUS "TTK_BUILD_STANDALONE_APPS: ${TTK_BUILD_STANDALONE_APPS}")
	message(STATUS "TTK_BUILD_VTK_WRAPPERS: ${TTK_BUILD_VTK_WRAPPERS}")
	if(TTK_BUILD_VTK_WRAPPERS)
		message(STATUS "  VTK_DIR: ${VTK_DIR}")
	endif()
	message(STATUS "ttk install -----------------------------------------------------------------")
	message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
	message(STATUS "TTK_INSTALL_PLUGIN_DIR: ${TTK_INSTALL_PLUGIN_DIR}")
	message(STATUS "-----------------------------------------------------------------------------")

endfunction()

option(TTK_BUILD_DOCUMENTATION "Build doxygen developer documentation" OFF)
if(TTK_BUILD_DOCUMENTATION)
	find_package(Doxygen)
	if (DOXYGEN_FOUND)
		configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ttk.doxygen
			${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen)
		add_custom_target(doc ALL
			${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			COMMENT "Generating API documentation with Doxygen" VERBATIM)
		install(DIRECTORY
			${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION
			${CMAKE_INSTALL_PREFIX}/share/doc/ttk)
		install(DIRECTORY
          ${CMAKE_SOURCE_DIR}/doc/img DESTINATION
			${CMAKE_INSTALL_PREFIX}/share/doc/ttk)
	endif()
endif()


find_package(ZLIB)
if(NOT ZLIB_INCLUDE_DIR)
  option(TTK_ENABLE_ZLIB "Enable Zlib support" OFF)
  message(WARNING "Zlib not found, disabling Zlib support in TTK.")
else()
  option(TTK_ENABLE_ZLIB "Enable Zlib support" ON)
endif()

# START_FINDZFP
find_path(ZFP_INCLUDE_DIR zfp.h
    PATHS
    $ENV{ZFP_ROOT_DIR}/include
    /opt/zfp/include)
set(ZFP_NAMES ${ZFP_NAMES} zfp)
find_library(ZFP_LIBRARY
    NAMES ${ZFP_NAMES}
    PATHS $ENV{ZFP_ROOT_DIR}/lib /opt/zfp/lib)
if (ZFP_LIBRARY AND ZFP_INCLUDE_DIR)
  set(ZFP_LIBRARIES ${ZFP_LIBRARY})
  set(ZFP_FOUND "YES")
else (ZFP_LIBRARY AND ZFP_INCLUDE_DIR)
  set(ZFP_FOUND "NO")
endif (ZFP_LIBRARY AND ZFP_INCLUDE_DIR)
if (ZFP_FOUND)
  if (NOT ZFP_FIND_QUIETLY)
    message(STATUS "Found ZFP: ${ZFP_LIBRARIES}")
  endif (NOT ZFP_FIND_QUIETLY)
else (ZFP_FOUND)
  if (ZFP_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find ZFP library...")
  endif (ZFP_FIND_REQUIRED)
endif (ZFP_FOUND)
mark_as_advanced(
    ZFP_LIBRARY
    ZFP_INCLUDE_DIR)
# END_FINDZFP
if (ZFP_FOUND)
  option(TTK_ENABLE_ZFP "Enable ZFP support" ON)
else()
  option(TTK_ENABLE_ZFP "Enable ZFP support" OFF)
  message(WARNING "ZFP not found, disabling ZFP support in TTK!")
endif()

if (NOT APPLE)
	if (MSVC)
  		option(TTK_ENABLE_OPENMP "Enable OpenMP support" FALSE)
	else()
  		option(TTK_ENABLE_OPENMP "Enable OpenMP support" TRUE)
	endif()
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
