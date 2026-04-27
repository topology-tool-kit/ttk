ttk_add_base_library(postProcessingTracking
  SOURCES
    postProcessingTracking.cpp
  HEADERS
    postProcessingTracking.h
  DEPENDS
    triangulation
    localizedTopologicalSimplification
    exTreeM
    ftmTreePP
    pathCompression
)

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(postProcessingTracking PUBLIC TTK_ENABLE_EIGEN)
  target_link_libraries(postProcessingTracking PUBLIC Eigen3::Eigen)
endif()
