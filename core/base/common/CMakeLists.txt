ttk_add_base_library(common
    SOURCES
        BaseClass.cpp
        Debug.cpp
        Os.cpp
    HEADERS
        BaseClass.h
        Cache.h
        CommandLineParser.h
        Debug.h
        DataTypes.h
        FlatJaggedArray.h
        MPIUtils.h
        OpenMP.h
        OrderDisambiguation.h
        Os.h
        ProgramBase.h
        Shuffle.h
        Timer.h
        VisitedMask.h
        Wrapper.h
        )

if(TTK_REDUCE_TEMPLATE_INSTANTIATIONS)
  # in the base layer, only this module is affected
  target_compile_definitions(common PRIVATE TTK_REDUCE_TEMPLATE_INSTANTIATIONS)
endif()
