ttk_add_baseCode_package(abstractTriangulation)
ttk_add_baseCode_package(zeroSkeleton)
ttk_add_baseCode_package(oneSkeleton)
ttk_add_baseCode_package(twoSkeleton)
ttk_add_baseCode_package(threeSkeleton)

# if the package is not a template, uncomment the following line
ttk_wrapup_library(libExplicitTriangulation "ExplicitTriangulation.cpp")

