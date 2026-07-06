ttk_add_base_library(timeTrackingPostProcessing
  SOURCES
    TimeTrackingPostProcessing.cpp
  HEADERS
    TimeTrackingPostProcessing.h
  DEPENDS
	  exTreeM
    ftmTreePP
    geometry
    localizedTopologicalSimplification
    pathCompression
    triangulation
)

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(timeTrackingPostProcessing PUBLIC TTK_ENABLE_EIGEN)
  target_link_libraries(timeTrackingPostProcessing PUBLIC Eigen3::Eigen)
endif()
