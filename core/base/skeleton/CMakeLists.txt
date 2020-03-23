if ("${ParaView_VERSION}" VERSION_GREATER_EQUAL "5.8.0" OR "${VTK_VERSION}" VERSION_GREATER_EQUAL "9.0.0")
  set(CellArrayHeader "CellArray_new.h")
  set(CellArrayDefine "CELL_ARRAY_NEW")
else()
  set(CellArrayHeader "CellArray_legacy.h")
  set(CellArrayDefine "CELL_ARRAY_LEGACY")
endif()
configure_file(${CellArrayHeader} CellArray.h COPYONLY)

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
  LINK
    common
  )

target_include_directories(skeleton PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_compile_definitions(skeleton PUBLIC ${CellArrayDefine})
