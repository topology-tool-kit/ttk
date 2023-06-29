ttk_add_base_library(constrainedGradientDescent
  SOURCES
    ConstrainedGradientDescent.cpp
  HEADERS
    ConstrainedGradientDescent.h
  DEPENDS
    common
    persistenceDiagram
)

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(constrainedGradientDescent PRIVATE TTK_ENABLE_EIGEN)
  target_link_libraries(constrainedGradientDescent PRIVATE Eigen3::Eigen)
endif()