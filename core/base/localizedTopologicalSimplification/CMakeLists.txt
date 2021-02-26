ttk_add_base_template_library(localizedTopologicalSimplification
  HEADERS
    LocalizedTopologicalSimplification.h
    SSCPropagation.h
  DEPENDS
    triangulation
    Boost::boost
)