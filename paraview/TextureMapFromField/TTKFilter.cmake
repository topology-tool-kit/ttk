option(TTK_BUILD_TEXTUREMAPFROMFIELD_FILTER "Build the TextureMapFromField filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_TEXTUREMAPFROMFIELD_FILTER)

if(${TTK_BUILD_TEXTUREMAPFROMFIELD_FILTER})
  ttk_register_pv_filter(ttkTextureMapFromField TextureMapFromField.xml)
endif()
