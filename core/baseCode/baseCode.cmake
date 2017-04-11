# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

  if(NOT BASECODE)

  # avoid some issues under windows
  if(COMMAND cmake_policy)
    cmake_policy(SET CMP0002 OLD)
    cmake_policy(SET CMP0003 NEW)
    cmake_policy(SET CMP0004 OLD)
    cmake_policy(SET CMP0011 NEW)
    cmake_policy(SET CMP0012 NEW)
    cmake_policy(SET CMP0023 OLD)
    cmake_policy(SET CMP0028 OLD)
    cmake_policy(SET CMP0054 NEW)
  endif(COMMAND cmake_policy)
  
  # kamikaze code compilation mode
  option(withKamikaze "Enable Kamikaze compilation mode" false)

  if(withKamikaze)
    set(PROJECT_FLAGS "${PROJECT_FLAGS} -DwithKamikaze"
      CACHE INTERNAL "PROJECT_FLAGS")
  endif(withKamikaze)

  # cross build
  # needs multilib support on the system (including vtk and paraview)
  # NOTE: gentoo does not provide multilib support for vtk or paraview
  # ...for now. let's be patient
  #  option(with32bits "Enable 32-bit mode" false)
  #if(with32bits)
  #  set(PROJECT_FLAGS "${PROJECT_FLAGS} -m32")
  #  set(PROJECT_DEP "${PROJECT_DEP} -m32")
  #endif(with32bits)

  # cpuOptimization - to disable for release on other computers
  option(withCpuOptimization "Enable CPU optimization (to disable for release)"
    true)
  if(withCpuOptimization)
    set(PROJECT_FLAGS "${PROJECT_FLAGS} -march=native")
  endif(withCpuOptimization)

  # openMP compilation flags
  option(withOpenMP "Enable OpenMP support" true)
  if(APPLE)
    # disable OpenMP by default under MacOS
    # TODO JAL Make this better
    set(withOpenMP false)
  endif(APPLE)
  if(withOpenMP)
    find_package(OpenMP )
    if(OPENMP_FOUND)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(PROJECT_FLAGS "${PROJECT_FLAGS} -DwithOpenMP"
        CACHE INTERNAL "PROJECT_FLAGS")
    endif(OPENMP_FOUND)
  endif(withOpenMP)


  # MPI compilation flags
  option(withMPI "Enable MPI support" false)
  if (withMPI)
    find_package(MPI )
    include_directories(${MPI_INCLUDE_PATH})
    set(PROJECT_FLAGS "${PROJECT_FLAGS} -DwithMPI"
      CACHE INTERNAL "PROJECT_FLAGS")
  endif(withMPI)

  message(STATUS
"baseCode --------------------------------------------------------------------")
#  message(STATUS
#    "32 bit compilation mode: ${with32bits}   (-Dwith32bits=)")
  message(STATUS
    "CPU optimizations: ${withCpuOptimization}   (-DwithCpuOptimization=)")
  message(STATUS 
    "Kamikaze compilation mode: ${withKamikaze}   (-DwithKamikaze=)")
  message(STATUS "OpenMP: ${withOpenMP}   (-DwithOpenMP=)")
  message(STATUS "MPI: ${withMPI}   (-DwithMPI=)")
  message(STATUS
"-----------------------------------------------------------------------------")

  # environment specific settings
  if(MSVC)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
#    set_target_properties(${PROJECT}
#      PROPERTIES
#      RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR})
#    set_target_properties(${PROJECT}
#      PROPERTIES
#      RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR})
#    set(PROJECT_FLAGS ${PROJECT_FLAGS} "/MP" CACHE INTERNAL "PROJECT_FLAGS")
  else(MSVC)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    if(COMPILER_SUPPORTS_CXX11)
      set(CXX_FLAGS "-std=c++11")
    else(COMPILER_SUPPORTS_CXX11)
      set(CXX_FLAGS "-std=c++0x")
    endif(COMPILER_SUPPORTS_CXX11)

    set(PROJECT_FLAGS
      "${PROJECT_FLAGS} -Wall -fPIC ${CXX_FLAGS}"
      CACHE INTERNAL "PROJECT_FLAGS")
    if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
      set (PROJECT_FLAGS
        "${PROJECT_FLAGS} -O3"
        CACHE INTERNAL "PROJECT_FLAGS")
    endif(NOT CMAKE_BUILD_TYPE MATCHES Debug)
  endif()

endif(NOT BASECODE)
