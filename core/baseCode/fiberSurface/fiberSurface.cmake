ttk_add_baseCode_package(geometry)
ttk_add_baseCode_package(triangulation)
ttk_add_optional_baseCode_package(rangeDrivenOctree)

# if the package is not a template, uncomment the following line
ttk_wrapup_library(libFiberSurface "FiberSurface.cpp")
