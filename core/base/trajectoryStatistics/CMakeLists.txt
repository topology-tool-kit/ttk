ttk_add_base_library(trajectoryStatistics
  SOURCES
    TrajectoryStatistics.cpp
  HEADERS
    TrajectoryStatistics.h
  DEPENDS
    triangulation
)

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(trajectoryStatistics PUBLIC TTK_ENABLE_EIGEN)
  target_link_libraries(trajectoryStatistics PUBLIC Eigen3::Eigen)
endif()
