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
        ArrayLinkedList.h   
        )

if(TTK_REDUCE_TEMPLATE_INSTANTIATIONS)
  # in the base layer, only this module is affected
  target_compile_definitions(common PRIVATE TTK_REDUCE_TEMPLATE_INSTANTIATIONS)
endif()

if(CMAKE_SIZEOF_VOID_P LESS_EQUAL 4)
  # provides compatibility between LongSimplexId and vtkIdType on i386
  # c.f. VTK_USE_64BIT_IDS in $PARAVIEW_SOURCE_DIR/VTK/Common/Core/CMakeLists.txt
  target_compile_definitions(common PUBLIC TTK_HW_IS_32BITS)
endif()

if(WIN32)
  # disable Windows warnings on some POSIX calls being insecure
  # (Windows provides secure but non-portable alternatives)
  target_compile_definitions(common PUBLIC _CRT_SECURE_NO_WARNINGS)
endif()
