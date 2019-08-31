option(TTK_PROGRAMBASE_FILTER "Build the ProgramBase filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_PROGRAMBASE_FILTER)

if(${TTK_PROGRAMBASE_FILTER})
  ttk_register_pv_module(ttkProgramBase)
endif()
