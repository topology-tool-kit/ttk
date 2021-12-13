# --- Prerequisites

set(CMAKE_CXX_STANDARD 11)

# --- Global Options

option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" ON)
cmake_dependent_option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" ON "TTK_BUILD_VTK_WRAPPERS" OFF)
option(TTK_BUILD_STANDALONE_APPS "Build the TTK Standalone Applications" ON)
option(TTK_WHITELIST_MODE "Explicitely enable each filter" OFF)
mark_as_advanced(TTK_WHITELIST_MODE BUILD_SHARED_LIBS)

# This option allows library to be built dynamic
# like the TopologyToolKit.so file for paraview
option(BUILD_SHARED_LIBS "Build TTK as shared lib" ON)

if(TTK_BUILD_STANDALONE_APPS AND NOT TTK_BUILD_VTK_WRAPPERS)
  message(WARNING "Can't build standalones without the VTK wrappers: disable")
  set(TTK_BUILD_STANDALONE_APPS OFF CACHE BOOL "Build the cmd and gui commands" FORCE)
endif()

# Set a default build type if none was specified
get_property(generator_is_multi_config GLOBAL
  PROPERTY GENERATOR_IS_MULTI_CONFIG)
if (NOT CMAKE_BUILD_TYPE AND NOT generator_is_multi_config)
  message(STATUS "Setting build type to 'Release'.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif ()

# the 'TTK_CELL_ARRAY_LAYOUT' drive the behavior of TTK to store cell points.
# This variable has two possible values "SingleArray" and "OffsetAndConnectiviy".
# * "SingleArray" use a layout compatible with VTK < 9 were the cell array store the
#   cells and their connectivity in a flat array
# * "OffsetAndConnectivity" use a layout comatible with VTK >= 9, having two arrays
#   (see https://vtk.org/doc/nightly/html/classvtkCellArray.html#details for more info)
set(TTK_CELL_ARRAY_LAYOUT "SingleArray" CACHE STRING "Layout for the cell array.")
set_property(CACHE TTK_CELL_ARRAY_LAYOUT PROPERTY STRINGS "SingleArray" "OffsetAndConnectivity")
mark_as_advanced(TTK_CELL_ARRAY_LAYOUT)

# issue #605
# workaround https://gitlab.kitware.com/paraview/paraview/-/issues/20324
option(TTK_ENABLE_MPI "Enable MPI support" FALSE)
if (TTK_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

if(TTK_BUILD_PARAVIEW_PLUGINS OR TTK_BUILD_VTK_WRAPPERS)
  # Find ParaView, otherwise VTK
  find_package(ParaView)
  if(ParaView_FOUND)
    # hande version manually so we do not have to include
    # files from VTK / ParaView in the code.
    # this is necessary to work with build folder directly (MacOS)
    add_definitions(-DPARAVIEW_VERSION_MAJOR=${ParaView_VERSION_MAJOR})
    add_definitions(-DPARAVIEW_VERSION_MINOR=${ParaView_VERSION_MINOR})
    # Layout to use for the CellArray is driven by ParaView version:
    # TODO: we can even hide the option here as the user should not change it in this case.
    if ("${ParaView_VERSION}" VERSION_GREATER_EQUAL "5.8.0")
      set(TTK_CELL_ARRAY_LAYOUT "OffsetAndConnectivity" CACHE STRING "Layout for the cell array." FORCE)
    else()
      set(TTK_CELL_ARRAY_LAYOUT "SingleArray" CACHE STRING "Layout for the cell array." FORCE)
    endif()
  else()
    find_package(VTK)
    if(VTK_FOUND)
      # hande version manually so we do not have to include
      # files from VTK / ParaView in the code.
      # this is necessary to work with build folder directly (MacOS)
      add_definitions(-DVTK_VERSION_MAJOR=${VTK_VERSION_MAJOR})
      add_definitions(-DVTK_VERSION_MINOR=${VTK_VERSION_MINOR})

      # Layout to use for the CellArray is driven by VTK version:
      # TODO: we can even hide the option here as the user should not change it in this case.
      if ("${VTK_VERSION}" VERSION_GREATER_EQUAL "9.0")
        set(TTK_CELL_ARRAY_LAYOUT "OffsetAndConnectivity" CACHE STRING "Layout for the cell array." FORCE)
      else()
        # VTK 8.2 is not supported, 8.90 is not official. Troubles incoming if you try this.
        # this should only works with VTK version close to the PV 5.7 internal one:
        # tag e4e8a4df9cc67fd2bb3dbb3b1c50a25177cbfe68
        message(WARNING "This VTK version is not supported. You may have compilation error.")
        set(TTK_CELL_ARRAY_LAYOUT "SingleArray" CACHE STRING "Layout for the cell array." FORCE)
      endif()
    endif()
  endif()
endif()

if(TTK_BUILD_PARAVIEW_PLUGINS)
  if(NOT ParaView_FOUND)
    message(FATAL_ERROR "TTK_BUILD_PARAVIEW_PLUGINS requires ParaView.")
  endif()
elseif(TTK_BUILD_VTK_WRAPPERS)
  if(NOT ParaView_FOUND AND NOT VTK_FOUND)
    message(FATAL_ERROR "TTK_BUILD_VTK_WRAPPERS requires ParaView or VTK.")
  endif()
endif()

option(TTK_ENABLE_64BIT_IDS "Enable processing on large datasets" OFF)
mark_as_advanced(TTK_ENABLE_64BIT_IDS)

option(TTK_ENABLE_KAMIKAZE "Enable Kamikaze compilation mode" OFF)
mark_as_advanced(TTK_ENABLE_KAMIKAZE)

option(TTK_ENABLE_CPU_OPTIMIZATION "Enable native CPU optimizations" ON)
mark_as_advanced(TTK_ENABLE_CPU_OPTIMIZATION)

option(TTK_ENABLE_DOUBLE_TEMPLATING "Use double templating for bivariate data" OFF)
mark_as_advanced(TTK_ENABLE_DOUBLE_TEMPLATING)

option(TTK_ENABLE_SHARED_BASE_LIBRARIES "Generate shared base libraries instead of static ones" ON)
mark_as_advanced(TTK_ENABLE_SHARED_BASE_LIBRARIES)
if(TTK_ENABLE_SHARED_BASE_LIBRARIES AND MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

option(TTK_BUILD_DOCUMENTATION "Build doxygen developer documentation" OFF)
if(TTK_BUILD_DOCUMENTATION)
  find_package(Doxygen)
endif()

# --- Dependencies

# for additional find packages
list(INSERT CMAKE_MODULE_PATH 0
  "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

# mandatory packages

find_package(Boost REQUIRED)
if(Boost_FOUND)
  message(STATUS "Found Boost ${Boost_VERSION} (${Boost_INCLUDE_DIR})")
endif()

# optional pakages

find_package(ZLIB QUIET)
if(ZLIB_FOUND)
  option(TTK_ENABLE_ZLIB "Enable Zlib support" ON)
  message(STATUS "Found Zlib ${ZLIB_VERSION_STRING} (${ZLIB_LIBRARIES})")
else()
  option(TTK_ENABLE_ZLIB "Enable Zlib support" OFF)
  message(STATUS "Zlib not found, disabling Zlib support in TTK.")
endif()

find_package(EMBREE 3.4 QUIET)
if(EMBREE_FOUND)
  option(TTK_ENABLE_EMBREE "Enable embree raytracing for ttkCinemaImaging" ON)
  message(STATUS "Found Embree ${EMBREE_VERSION} (${EMBREE_LIBRARY})")
else()
  option(TTK_ENABLE_EMBREE "Enable embree raytracing for ttkCinemaImaging" OFF)
  message(STATUS "EMBREE library not found, disabling embree support in TTK.")
endif()

find_package(Graphviz QUIET)
if(Graphviz_FOUND)
  option(TTK_ENABLE_GRAPHVIZ "Enable GraphViz support" ON)
  message(STATUS "Found GraphViz (${Graphviz_CGRAPH_LIBRARY})")
else()
  option(TTK_ENABLE_GRAPHVIZ "Enable GraphViz support" OFF)
  message(STATUS "GraphViz not found, disabling GraphViz support in TTK.")
endif()

find_package(SQLite3 QUIET)
if(SQLite3_FOUND)
  option(TTK_ENABLE_SQLITE3 "Enable SQLITE3 support" ON)
  message(STATUS "Found SQLite3 ${SQLite3_VERSION} (${SQLite3_LIBRARIES})")
else()
  option(TTK_ENABLE_SQLITE3 "Enable SQLITE3 support" OFF)
  message(STATUS "SQLite3 not found, disabling SQLite3 support in TTK.")
endif()

find_package(ZFP QUIET)
if(ZFP_INCLUDE_DIRS)
  option(TTK_ENABLE_ZFP "Enable ZFP support" ON)
  get_property(ZFP_LIB TARGET ${ZFP_LIBRARIES} PROPERTY IMPORTED_LOCATION_RELEASE)
  message(STATUS "Found ZFP ${ZFP_VERSION} (${ZFP_LIB})")
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

find_package(Eigen3 3.3 QUIET NO_MODULE)
if(EIGEN3_FOUND)
  option(TTK_ENABLE_EIGEN "Enable Eigen3 support" ON)
  message(STATUS "Found Eigen ${Eigen3_VERSION} (${EIGEN3_INCLUDE_DIR})")

  find_package(Spectra 1.0.0 QUIET)
  if(Spectra_FOUND)
    option(TTK_ENABLE_SPECTRA "Enable Spectra support" ON)
    get_property(SPECTRA_INCLUDES TARGET Spectra::Spectra PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "Found Spectra ${Spectra_VERSION} (${SPECTRA_INCLUDES})")
  else()
    option(TTK_ENABLE_SPECTRA "Enable Spectra support" OFF)
    message(STATUS "Spectra >=1.0.0 not found, disabling Spectra support in TTK.")
  endif()
else()
  option(TTK_ENABLE_EIGEN "Enable Eigen3 support" OFF)
  message(STATUS "Eigen not found, disabling Eigen support in TTK.")

  option(TTK_ENABLE_SPECTRA "Enable Spectra support" OFF)
  message(STATUS "Spectra not found, disabling Spectra support in TTK.")
endif()

find_package(Python3 COMPONENTS Development NumPy QUIET)
if(Python3_FOUND AND Python3_NumPy_FOUND)
  include_directories(SYSTEM ${Python3_INCLUDE_DIRS})
  set(TTK_PYTHON_MAJOR_VERSION "${Python3_VERSION_MAJOR}"
    CACHE INTERNAL "TTK_PYTHON_MAJOR_VERSION")
  set(TTK_PYTHON_MINOR_VERSION "${Python3_VERSION_MINOR}"
    CACHE INTERNAL "TTK_PYTHON_MINOR_VERSION")

  option(TTK_ENABLE_SCIKIT_LEARN "Enable scikit-learn support" ON)
  message(STATUS "Found Python ${Python3_VERSION} (${Python3_EXECUTABLE})")
else()
  option(TTK_ENABLE_SCIKIT_LEARN "Enable scikit-learn support" OFF)
  message(STATUS "Improper Python/NumPy setup. Disabling scikit-learn support in TTK.")
endif()

if(MSVC)
  option(TTK_ENABLE_OPENMP "Enable OpenMP support" FALSE)
else()
  option(TTK_ENABLE_OPENMP "Enable OpenMP support" TRUE)
endif()
if(TTK_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
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


find_package(WEBSOCKETPP QUIET)
if(WEBSOCKETPP_FOUND)
  option(TTK_ENABLE_WEBSOCKETPP "Enable WebSocketIO module" ON)
  message(STATUS "Found WebSocketPP ${WEBSOCKETPP_VERSION} (${WEBSOCKETPP_INCLUDE_DIR}), enabling WebSocketIO module in TTK.")
else()
  option(TTK_ENABLE_WEBSOCKETPP "Enable WebSocketIO module" OFF)
  message(STATUS "WebSocketPP not found, disabling WebSocketIO module in TTK.")
endif()

# --- Install path

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

# ParaView plugins go to a subdirectory with this name
set(TTK_PLUGIN_SUBDIR "TopologyToolKit")

# Install rpath
if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH TRUE)
endif()
if(NOT DEFINED CMAKE_BUILD_WITH_INSTALL_RPATH)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()
if(NOT DEFINED CMAKE_INSTALL_RPATH)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()
if(NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()
