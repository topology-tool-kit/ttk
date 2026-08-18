ttk_add_base_library(trackingPostProcessing
  SOURCES
    TrackingPostProcessing.cpp
  HEADERS
    TrackingPostProcessing.h
  DEPENDS
	  exTreeM
    ftmTreePP
    geometry
    localizedTopologicalSimplification
    pathCompression
    triangulation
)

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(trackingPostProcessing PUBLIC TTK_ENABLE_EIGEN)
  target_link_libraries(trackingPostProcessing PUBLIC Eigen3::Eigen)
endif()
