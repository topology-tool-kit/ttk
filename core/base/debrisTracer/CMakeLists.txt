ttk_add_base_library(debrisTracer
  SOURCES
    DebrisTracer.cpp
  HEADERS
    DebrisTracer.h
  DEPENDS
    triangulation
	persistenceDiagram
	topologicalSimplification
	exTreeM
	pathCompression
)

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(debrisTracer PUBLIC TTK_ENABLE_EIGEN)
  target_link_libraries(debrisTracer PUBLIC Eigen3::Eigen)
endif()

