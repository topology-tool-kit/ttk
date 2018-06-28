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

	set(CMAKE_SKIP_BUILD_RPATH TRUE)
	set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
	set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
	set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

	add_library(${library} STATIC ${ARG_SOURCES})
	set_target_properties(${library} PROPERTIES
		POSITION_INDEPENDENT_CODE TRUE)

	target_include_directories(${library} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:include/ttk/base>)

	if (ARG_LINK)
		target_link_libraries(${library} PUBLIC ${ARG_LINK})
	endif()

	ttk_set_compile_options(${library})

	install(TARGETS ${library}
		EXPORT TTKBaseConfig
		RUNTIME DESTINATION bin/ttk
		ARCHIVE DESTINATION lib/ttk
		LIBRARY DESTINATION lib/ttk)

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

	if (TTK_ENABLE_KAMIKAZE)
		target_compile_definitions(${library} INTERFACE TTK_ENABLE_KAMIKAZE)
	endif()

	string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)
  if (uppercase_CMAKE_BUILD_TYPE MATCHES RELEASE)
    if (TTK_ENABLE_CPU_OPTIMIZATION)
      target_compile_options(${library} INTERFACE -march=native -O3)
    endif()
  endif()

	if (TTK_ENABLE_OPENMP)
		target_compile_definitions(${library} INTERFACE TTK_ENABLE_OPENMP)
		target_compile_options(${library} INTERFACE ${OpenMP_CXX_FLAGS})
	endif()

	if (ARG_LINK)
		target_link_libraries(${library} INTERFACE ${ARG_LINK})
	endif()

	if (TTK_ENABLE_MPI)
		target_link_libraries(${library} INTERFACE ${MPI_CXX_LIBRARIES})
	endif()

  if (TTK_ENABLE_ZFP)
		target_compile_definitions(${library} INTERFACE TTK_ENABLE_ZFP)
		target_include_directories(${library} INTERFACE ${ZFP_INCLUDE_DIR})
	endif()
	
	if (TTK_ENABLE_ZLIB)
		target_compile_definitions(${library} INTERFACE TTK_ENABLE_ZLIB)
		target_include_directories(${library} INTERFACE ${ZLIB_INCLUDE_DIR})
	endif()

	
	install(TARGETS ${library}
		EXPORT TTKBaseConfig
		RUNTIME DESTINATION bin/ttk
		ARCHIVE DESTINATION lib/ttk
		LIBRARY DESTINATION lib/ttk)

	install(FILES ${ARG_SOURCES} DESTINATION include/ttk/base)
endfunction()

set(LIB_LIST "")
file(GLOB BASECODE_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)
foreach(BASECODE ${BASECODE_DIRS})
	if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${BASECODE})
		add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${BASECODE})
		set(LIB_LIST "${BASECODE};${LIB_LIST}")
	endif()
endforeach()

add_library(baseAll INTERFACE)
target_link_libraries(baseAll INTERFACE ${LIB_LIST})

install(TARGETS baseAll
	EXPORT TTKBaseConfig
	RUNTIME DESTINATION bin/ttk
	ARCHIVE DESTINATION lib/ttk
	LIBRARY DESTINATION lib/ttk)

install(EXPORT TTKBaseConfig
	NAMESPACE ttk::base::
	DESTINATION lib/cmake/ttk)

