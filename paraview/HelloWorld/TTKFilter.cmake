# Allows to disable each filter
option(TTK_BUILD_HELLOWORLD_FILTER "Build the HelloWorld filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_HELLOWORLD_FILTER)

if(${TTK_BUILD_HELLOWORLD_FILTER})
  ttk_register_pv_filter(ttkHelloWorld HelloWorld.xml)
endif()
