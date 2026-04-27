ttk_add_base_library(postProcessingTracking
  SOURCES
    PostProcessingTracking.cpp
  HEADERS
    PostProcessingTracking.h
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
