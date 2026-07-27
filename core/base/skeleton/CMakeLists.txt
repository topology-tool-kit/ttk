if ("${TTK_CELL_ARRAY_LAYOUT}" STREQUAL "OffsetAndConnectivity")
  set(CellArrayHeader "CellArray_new.h")
  set(CellArrayDefine "TTK_CELL_ARRAY_NEW")
elseif("${TTK_CELL_ARRAY_LAYOUT}" STREQUAL "SingleArray")
  set(CellArrayHeader "CellArray_legacy.h")
  set(CellArrayDefine "TTK_CELL_ARRAY_LEGACY")
else()
  message(FATAL_ERROR "Invalid TTK_CELL_ARRAY_LAYOUT.")
endif()

configure_file(${CellArrayHeader} CellArray.h COPYONLY)
install(
  FILES
   "${CMAKE_CURRENT_BINARY_DIR}/CellArray.h"
  DESTINATION
    "include/ttk/base"
  )

ttk_add_base_library(skeleton
  SOURCES
    CellArray.cpp
    ZeroSkeleton.cpp
    OneSkeleton.cpp
    TwoSkeleton.cpp
    ThreeSkeleton.cpp
  HEADERS
    ${CellArrayHeader}
    ZeroSkeleton.h
    OneSkeleton.h
    TwoSkeleton.h
    ThreeSkeleton.h
  DEPENDS
    common
  )

target_include_directories(skeleton PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_compile_definitions(skeleton PUBLIC ${CellArrayDefine})
target_link_libraries(skeleton PRIVATE Boost::boost)
