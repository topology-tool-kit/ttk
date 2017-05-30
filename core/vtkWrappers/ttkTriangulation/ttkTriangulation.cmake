ttk_add_baseCode_package(triangulation)

# This class is not meant to be used by a paraview plugin, use library.
ttk_wrapup_library(libttkTriangulation "ttkTriangulation.cpp")