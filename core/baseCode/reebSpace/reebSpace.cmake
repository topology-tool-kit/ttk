ttk_add_baseCode_package(fiberSurface)
ttk_add_baseCode_package(geometry)
ttk_add_baseCode_package(jacobiSet)
ttk_add_baseCode_package(triangulation)

# if the package is not a template, uncomment the following line
ttk_wrapup_library(libReebSpace "ReebSpace.cpp")
