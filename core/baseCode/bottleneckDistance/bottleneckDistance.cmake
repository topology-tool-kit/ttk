ttk_add_baseCode_package(triangulation)
ttk_add_baseCode_package(persistenceDiagram)

# if the package is a pure template class, comment the following line
ttk_wrapup_library(libBottleneckDistance "BottleneckDistance.cpp Munkres.cpp")
