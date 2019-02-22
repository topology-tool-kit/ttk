ttk_add_base_library(harmonicFieldComputation
                     SOURCES
                     HarmonicFieldComputation.cpp
                     HEADERS
                     HarmonicFieldComputation.h
                     LINK
                     geometry
                     triangulation)

target_compile_options(harmonicFieldComputation PRIVATE -Wextra -pedantic)
