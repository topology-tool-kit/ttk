ttk_add_base_library(common
    SOURCES
        BaseClass.cpp
        BaseMPIClass.cpp
        Debug.cpp
        Os.cpp
    HEADERS
        BaseClass.h
        BaseMPIClass.h
        CommandLineParser.h
        Debug.h
        DataTypes.h
        FlatJaggedArray.h
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
