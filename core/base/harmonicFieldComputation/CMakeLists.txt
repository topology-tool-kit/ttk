ttk_add_base_library(harmonicFieldComputation
                     SOURCES
                     HarmonicFieldComputation.cpp
                     HEADERS
                     HarmonicFieldComputation.h
                     LINK
                     triangulation)

target_compile_options(harmonicFieldComputation PRIVATE -Wextra -pedantic)
