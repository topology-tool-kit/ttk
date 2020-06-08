set(CMAKE_CXX_STANDARD 11)

# Set a predefined build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release'.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# Options & dependencies

# CellArray use legacy single array by default
set(TTK_CELL_ARRAY_LAYOUT "SingleArray" CACHE STRING "Layout for the cell array." FORCE)
set_property(CACHE TTK_CELL_ARRAY_LAYOUT PROPERTY STRINGS "SingleArray" "OffsetAndConnectivity")
mark_as_advanced(TTK_CELL_ARRAY_LAYOUT)

# Find Paraview OR VTK if needed
if(TTK_BUILD_PARAVIEW_PLUGINS)
  find_package(ParaView REQUIRED)
  # hande version manually so we do not have to include
  # files from VTK / ParaView in the code.
  # this is necessary to work with build folder directly (MacOS)
  add_definitions(-DPARAVIEW_VERSION_MAJOR=${ParaView_VERSION_MAJOR})
  add_definitions(-DPARAVIEW_VERSION_MINOR=${ParaView_VERSION_MINOR})
  # Layout to use for the CellArray is driven by ParaView version:
  if ("${ParaView_VERSION}" VERSION_GREATER_EQUAL "5.8.0")
    set(TTK_CELL_ARRAY_LAYOUT "OffsetAndConnectivity" CACHE STRING "Layout for the cell array." FORCE)
  else()
    set(TTK_CELL_ARRAY_LAYOUT "SingleArray" CACHE STRING "Layout for the cell array." FORCE)
  endif()
  # TODO: we can even hide the option here as the user should not change it in this case.

elseif(TTK_BUILD_VTK_WRAPPERS)
  find_package(VTK REQUIRED)
  add_definitions(-DVTK_VERSION_MAJOR=${VTK_VERSION_MAJOR})
  add_definitions(-DVTK_VERSION_MINOR=${VTK_VERSION_MINOR})

  # Layout to use for the CellArray is driven by VTK version:
  if ("${VTK_VERSION}" VERSION_GREATER_EQUAL "9.0.0")
    set(TTK_CELL_ARRAY_LAYOUT "OffsetAndConnectivity" CACHE STRING "Layout for the cell array." FORCE)
  else()
    set(TTK_CELL_ARRAY_LAYOUT "SingleArray" CACHE STRING "Layout for the cell array." FORCE)
  endif()
  # TODO: we can even hide the option here as the user should not change it in this case.
endif()

option(TTK_ENABLE_64BIT_IDS "Enable processing on large datasets" OFF)
mark_as_advanced(TTK_ENABLE_64BIT_IDS)

option(TTK_ENABLE_KAMIKAZE "Enable Kamikaze compilation mode" OFF)
mark_as_advanced(TTK_ENABLE_KAMIKAZE)

option(TTK_ENABLE_CPU_OPTIMIZATION "Enable native CPU optimizations" ON)
mark_as_advanced(TTK_ENABLE_CPU_OPTIMIZATION)

option(TTK_BUILD_DOCUMENTATION "Build doxygen developer documentation" OFF)
if(TTK_BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ttk.doxygen
      ${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen)
    add_custom_target(doc
      ALL
        ${DOXYGEN_EXECUTABLE}
        ${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen
      WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT
        "Generating API documentation with Doxygen"
      VERBATIM
      )
    install(
      DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/doc/html
      DESTINATION
        ${CMAKE_INSTALL_PREFIX}/share/doc/ttk
        )
    install(
      DIRECTORY
        ${CMAKE_SOURCE_DIR}/doc/img
      DESTINATION
        ${CMAKE_INSTALL_PREFIX}/share/doc/ttk
        )
  endif()
endif()

find_package(Boost COMPONENTS system)
if(NOT Boost_FOUND)
  find_package(Boost REQUIRED)
  if(Boost_FOUND)
    message(STATUS "BOOST_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
  endif()
endif()

find_package(ZLIB)
if(NOT ZLIB_FOUND)
  option(TTK_ENABLE_ZLIB "Enable Zlib support" OFF)
  message(STATUS "Zlib not found, disabling Zlib support in TTK.")
else()
  option(TTK_ENABLE_ZLIB "Enable Zlib support" ON)
endif()

# START_FIND_GRAPHVIZ
find_path(GRAPHVIZ_INCLUDE_DIR
  NAMES
    graphviz/cgraph.h
  HINTS
    ${_GRAPHVIZ_INCLUDE_DIR}
    )

find_library(GRAPHVIZ_CDT_LIBRARY
  NAMES
    cdt
  HINTS
    ${_GRAPHVIZ_LIBRARY_DIR}
    )

find_library(GRAPHVIZ_GVC_LIBRARY
  NAMES
    gvc
  HINTS
    ${_GRAPHVIZ_LIBRARY_DIR}
    )

find_library(GRAPHVIZ_CGRAPH_LIBRARY
  NAMES
    cgraph
  HINTS
    ${_GRAPHVIZ_LIBRARY_DIR}
    )

find_library(GRAPHVIZ_PATHPLAN_LIBRARY
  NAMES
    pathplan
  HINTS
    ${_GRAPHVIZ_LIBRARY_DIR}
    )

if(GRAPHVIZ_INCLUDE_DIR
    AND GRAPHVIZ_CDT_LIBRARY
    AND GRAPHVIZ_GVC_LIBRARY
    AND GRAPHVIZ_CGRAPH_LIBRARY
    AND GRAPHVIZ_PATHPLAN_LIBRARY
    )
  set(GRAPHVIZ_FOUND "YES")
else()
  set(GRAPHVIZ_FOUND "NO")
endif()

if(GRAPHVIZ_FOUND)
  option(TTK_ENABLE_GRAPHVIZ "Enable GraphViz support" ON)
  message(STATUS "GraphViz found")
else()
  option(TTK_ENABLE_GRAPHVIZ "Enable GraphViz support" OFF)
  message(STATUS "GraphViz not found, disabling GraphViz support in TTK.")
endif()
# END_FIND_GRAPHVIZ

# START_FIND_SQLITE3
find_path(SQLITE3_INCLUDE_DIR sqlite3.h
  PATHS
    /usr/include
    /usr/local/include
    )
set(SQLITE3_NAMES ${SQLITE3_NAMES} sqlite3)
find_library(SQLITE3_LIBRARY
  NAMES
    ${SQLITE3_NAMES}
  PATHS
    $ENV{SQLITE3_ROOT_DIR}/lib /opt/sqlite3/lib
    )
if (SQLITE3_LIBRARY AND SQLITE3_INCLUDE_DIR)
  set(SQLITE3_LIBRARIES ${SQLITE3_LIBRARY})
  set(SQLITE3_FOUND "YES")
else (SQLITE3_LIBRARY AND SQLITE3_INCLUDE_DIR)
  set(SQLITE3_FOUND "NO")
endif (SQLITE3_LIBRARY AND SQLITE3_INCLUDE_DIR)
if (SQLITE3_FOUND)
  if (NOT SQLITE3_FIND_QUIETLY)
    message(STATUS "Found SQLITE3: ${SQLITE3_LIBRARIES}")
  endif (NOT SQLITE3_FIND_QUIETLY)
else (SQLITE3_FOUND)
  if (SQLITE3_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find SQLITE3 library...")
  endif (SQLITE3_FIND_REQUIRED)
endif (SQLITE3_FOUND)
mark_as_advanced(SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR)
# END_FINDSQLITE3
if(SQLITE3_FOUND)
  option(TTK_ENABLE_SQLITE3 "Enable SQLITE3 support" ON)
else()
  option(TTK_ENABLE_SQLITE3 "Enable SQLITE3 support" OFF)
  message(STATUS "SQLITE3 not found, disabling SQLITE3 support in TTK.")
endif()

find_package(ZFP QUIET)
if(ZFP_INCLUDE_DIRS)
  option(TTK_ENABLE_ZFP "Enable ZFP support" ON)
else()
  option(TTK_ENABLE_ZFP "Enable ZFP support" OFF)
  message(STATUS "ZFP not found, disabling ZFP support in TTK.")
endif()
if(NOT TTK_ENABLE_ZFP)
  # we do not want ZFP targets to remains if ZFP disable.
  # a bit hacky but there is no clean way to remove the corresponding targets
  unset(ZFP_DIR CACHE)
  find_package(ZFP QUIET)
endif()

find_package(Eigen3 3.3 NO_MODULE)
if(EIGEN3_FOUND)
  option(TTK_ENABLE_EIGEN "Enable Eigen3 support" ON)

  find_path(SPECTRA_INCLUDE_DIR Spectra/SymEigsSolver.h)
  if(SPECTRA_INCLUDE_DIR STREQUAL "SPECTRA_INCLUDE_DIR-NOTFOUND")
    option(TTK_ENABLE_SPECTRA "Enable Spectra support" OFF)
    message(STATUS "Spectra not found, disabling Spectra support in TTK.")
  else()
    option(TTK_ENABLE_SPECTRA "Enable Spectra support" ON)
  endif()
else()
  option(TTK_ENABLE_EIGEN "Enable Eigen3 support" OFF)
  message(STATUS "Eigen not found, disabling Eigen support in TTK.")

  option(TTK_ENABLE_SPECTRA "Enable Spectra support" OFF)
  message(STATUS "Spectra not found, disabling Spectra support in TTK.")
endif()

# scikit-learn support is disabled by default for now under MacOs
if(APPLE)
  option(TTK_ENABLE_SCIKIT_LEARN "Enable scikit-learn support" OFF)
  message(STATUS "Disabling scikit-learn support by default under MacOs.")
endif()

if(NOT APPLE)
  if(MSVC)
    option(TTK_ENABLE_OPENMP "Enable OpenMP support" FALSE)
  else()
    option(TTK_ENABLE_OPENMP "Enable OpenMP support" TRUE)
  endif()
endif()
option(TTK_ENABLE_MPI "Enable MPI support" FALSE)

if(TTK_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OPENMP_FOUND)
    option(TTK_ENABLE_OMP_PRIORITY
      "Gives tasks priority, high perf improvement"
      OFF
      )
    if(OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 4
        AND OpenMP_CXX_VERSION_MINOR GREATER_EQUAL 5)
      set(TTK_ENABLE_OMP_PRIORITY
        OFF
        CACHE
        BOOL
        "Enable priorities on opnemp tasks"
        FORCE
        )
    endif()

    mark_as_advanced(TTK_ENABLE_OMP_PRIORITY)

  endif()
else()
  if(TTK_ENABLE_OMP_PRIORITY)
    # priorities are only meaningful when openmp is on
    set(TTK_ENABLE_OMP_PRIORITY
      OFF
      CACHE
      BOOL
      "Enable priorities on opnemp tasks"
      FORCE
      )
  endif()
endif()



if (TTK_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

# Install path

include(GNUInstallDirs)
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
endif()

# Install rapth

set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
